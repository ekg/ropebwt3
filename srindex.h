#ifndef RB3_SRINDEX_H
#define RB3_SRINDEX_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * SR-index: subsampled r-index with phi function and toehold tracking.
 *
 * The phi function maps SA[k] -> SA[k-1] and is piecewise linear over r
 * intervals (where r = number of BWT runs). We store the breakpoints in
 * sorted order for O(log r) binary search evaluation.
 *
 * The toehold is a known text position SA[hi] maintained during backward
 * search. It is updated when a BWT run boundary is crossed, using stored
 * SA samples at run boundaries.
 *
 * The subsampling parameter s controls a space-time tradeoff:
 *   s=1: full r-index, O(r) space, O(log r) per locate
 *   s>1: subsampled, O(r + n/s) space, O(s) per initial locate + O(log r) per phi step
 *
 * Reference:
 *   Gagie, Navarro, Prezza, JACM 2020 (r-index)
 *   Cobas et al., CPM 2021 (subsampled r-index)
 */

typedef struct {
	int64_t n_runs;       /* number of BWT runs (r) */
	int64_t n;            /* total BWT length */
	int32_t s;            /* subsampling parameter; s=1 means no subsampling */
	int64_t n_samples;    /* number of SA samples at run boundaries (up to 2r) */
	/*
	 * Phi function representation:
	 * phi_sa[0..n_runs-1]: SA values at the start of each BWT run, sorted by SA value.
	 * phi_da[0..n_runs-1]: For each sorted entry phi_sa[i], the SA value at position
	 *                      (BWT_pos - 1), i.e. the SA value of the preceding position.
	 *                      phi(phi_sa[i]) = phi_da[i] for the breakpoints.
	 *
	 * For a general SA value v in [phi_sa[i], phi_sa[i+1]):
	 *   phi(v) = phi_da[i] + (v - phi_sa[i])
	 * because phi is linear (SA[k]-1 = SA[k-1]) within a BWT run.
	 */
	int64_t *phi_sa;     /* sorted SA values at BWT run starts (breakpoints) */
	int64_t *phi_da;     /* phi values at breakpoints */
	/*
	 * Toehold support: SA samples at run boundaries indexed by BWT position.
	 * run_pos[i]:  BWT position of the last character in run i (0-indexed)
	 * run_sa[i]:   SA[run_pos[i]], the text position at end of run i
	 *
	 * During backward search, when hi crosses a run boundary, we look up
	 * the stored SA value from run_sa[].
	 */
	int64_t *run_pos;    /* BWT position of last char in each run */
	int64_t *run_sa;     /* SA value at each run_pos */
	/*
	 * Subsampled SA: BWT positions where SA[pos] % s == 0.
	 * For s=1, this contains all run boundary samples (= run_pos/run_sa copy).
	 * For s>1, this contains all BWT positions with SA divisible by s,
	 * collected during the backward walk from sentinels.
	 *
	 * sub_pos[] is sorted by BWT position for O(log(n/s)) binary search.
	 */
	int64_t n_sub;       /* number of subsampled SA entries */
	int64_t *sub_pos;    /* BWT positions, sorted */
	int64_t *sub_sa;     /* SA values at sub_pos positions */
	/*
	 * Multi-string support: cumulative sequence lengths and text-order mapping.
	 * For multi-string BWTs, SA values are absolute text positions.
	 * cum_len[i] = start of the i-th sequence in text order.
	 * cum_len[m] = total text length.
	 * text_order_sid[i] = BWT sentinel position (= sequence ID) for the i-th
	 *                     sequence in text order.
	 */
	int64_t m;           /* number of sentinels (sequences) */
	int64_t *cum_len;    /* cumulative lengths, size m+1 */
	int64_t *text_order_sid; /* sentinel BWT positions in text order, size m */
} rb3_srindex_t;

/* Build the SR-index phi function from an FM-index.
 * @param f          FM-index (must be loaded)
 * @param s          subsampling parameter (1 = no subsampling)
 * @param n_threads  number of threads for SA computation
 * @return           allocated SR-index, or NULL on failure
 */
rb3_srindex_t *rb3_srindex_build(const void *f, int32_t s, int n_threads);

/* Evaluate the phi function: phi(sa_val) = SA[k-1] where SA[k] = sa_val.
 * @param sr       SR-index
 * @param sa_val   a suffix array value (text position)
 * @return         SA[k-1], or -1 if sa_val is the minimum (SA[0])
 */
int64_t rb3_srindex_phi(const rb3_srindex_t *sr, int64_t sa_val);

/* Look up the toehold: given a BWT position that is at or near the end of
 * a run, return the stored SA value.
 * @param sr       SR-index
 * @param bwt_pos  BWT position (typically hi-1 from SA interval [lo,hi))
 * @return         SA[bwt_pos] if bwt_pos is at a run boundary, or -1
 */
int64_t rb3_srindex_toehold(const rb3_srindex_t *sr, int64_t bwt_pos);

/* Locate all occurrences in SA interval [lo, hi) given a toehold.
 * Uses phi function to enumerate positions starting from the toehold.
 * @param sr       SR-index
 * @param lo       start of SA interval
 * @param hi       end of SA interval
 * @param toehold_sa  SA[hi-1], a known text position (the toehold)
 * @param out      output array of size >= (hi - lo), filled with SA values
 * @return         number of positions written
 */
int64_t rb3_srindex_locate(const rb3_srindex_t *sr, int64_t lo, int64_t hi,
                           int64_t toehold_sa, int64_t *out);

/* Locate a single BWT position by walking LF until a subsampled SA sample.
 * @param sr       SR-index (must have s > 0)
 * @param f        FM-index (for LF-mapping)
 * @param bwt_pos  BWT position to locate
 * @return         SA[bwt_pos] (text position), or -1 on failure
 */
int64_t rb3_srindex_locate_one(const rb3_srindex_t *sr, const void *f, int64_t bwt_pos);

/* Locate all occurrences in SA interval [lo, hi).
 * Resolves the toehold automatically via locate_one, then uses phi.
 * @param sr       SR-index
 * @param f        FM-index (for LF-mapping in locate_one)
 * @param lo       start of SA interval
 * @param hi       end of SA interval
 * @param positions output array, filled with text positions
 * @param max_pos  maximum number of positions to return
 * @return         number of positions written, or -1 on failure
 */
int64_t rb3_srindex_locate_all(const rb3_srindex_t *sr, const void *f,
                               int64_t lo, int64_t hi, int64_t *positions, int64_t max_pos);

/* Free SR-index memory. */
void rb3_srindex_destroy(rb3_srindex_t *sr);

/* Serialize SR-index to a binary file (.sri format). */
int rb3_srindex_dump(const rb3_srindex_t *sr, const char *fn);

/* Deserialize SR-index from a binary file (.sri format). */
rb3_srindex_t *rb3_srindex_restore(const char *fn);

#ifdef __cplusplus
}
#endif

#endif
