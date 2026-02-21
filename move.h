#ifndef RB3_MOVE_H
#define RB3_MOVE_H

#include <stdint.h>
#include "fm-index.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rb3_move_row_s {
	int64_t p;       // Starting BWT offset of this run
	int64_t pi;      // LF[p] — LF-mapping of the run head
	int64_t xi;      // Index of the destination row containing pi
	int64_t len;     // Run length
	int16_t dist[RB3_ASIZE]; // Distance (in rows) to nearest run of each character; 0 = self
	int8_t c;        // BWT character of this run (0-5)
	uint8_t pad_[3]; // padding to 48 bytes for mmap alignment
} rb3_move_row_t;

typedef struct rb3_move_s {
	int64_t n_runs;          // Number of runs (r)
	int64_t bwt_len;         // Total BWT length (n)
	int64_t acc[7];          // Cumulative character counts C[]
	int32_t d;               // Run-splitting depth (0 = no splitting)
	rb3_move_row_t *rows;    // Flat array of rows
	int _fd;                 // file descriptor for mmap (internal)
	void *_mmap_base;        // mmap base address (internal, NULL if heap-allocated)
	size_t _mmap_len;        // mmap region length (internal)
} rb3_move_t;

// Build move table from FM-index (no splitting)
rb3_move_t *rb3_move_build(const rb3_fmi_t *f);
void rb3_move_destroy(rb3_move_t *m);

// Run splitting: split long runs so fast-forward is bounded by < 2d steps.
// Modifies m in place. d is the splitting depth (typically 2 or 3).
void rb3_move_split(rb3_move_t *m, int32_t d);

// Precompute reposition distances for each row.
void rb3_move_precompute_dist(rb3_move_t *m);

// LF-mapping: compute LF(pos) given current position and run index.
// Updates *run_idx to the destination run index.
int64_t rb3_move_lf(const rb3_move_t *m, int64_t pos, int64_t *run_idx);

// Reposition: jump to nearest run of character c from run_idx.
// Returns the new run index.
int64_t rb3_move_reposition(const rb3_move_t *m, int64_t run_idx, int8_t c);

// Combined reposition + LF-mapping: one backward search step.
// Given current pos in run run_idx, extend backward with character c.
// Updates *run_idx and returns the new BWT position.
int64_t rb3_move_step(const rb3_move_t *m, int64_t pos, int64_t *run_idx, int8_t c);

// Count occurrences of pattern[0..len-1] (encoded as 0-5 integers).
// Returns the SA interval size (number of occurrences).
int64_t rb3_move_count(const rb3_move_t *m, int len, const uint8_t *pattern);

struct rb3_lcp_s; /* forward declaration for LCP integration */

// Precompute per-move-row thresholds from LCP thresholds.
// Maps each move row to its corresponding LCP run threshold via linear scan.
// Returns malloc'd array of m->n_runs int64_t values. Caller must free().
int64_t *rb3_move_lcp_thresholds(const rb3_move_t *m, const struct rb3_lcp_s *lcp);

// Map each move row to its LCP run index.
// Returns malloc'd array of m->n_runs int64_t values. Caller must free().
int64_t *rb3_move_lcp_run_map(const rb3_move_t *m, const struct rb3_lcp_s *lcp);

// One backward step of matching statistics with move + LCP.
// pos: current BWT position in run *run_idx, *match_len: current match length.
// c: query character (nt6, 1-5). run_map: move-row-to-LCP-run mapping from
// rb3_move_lcp_run_map(). Returns new BWT position; updates *run_idx
// and *match_len. Returns -1 if c does not exist in the BWT.
int64_t rb3_move_ms_step(const rb3_move_t *m, const int64_t *run_map,
                         const struct rb3_lcp_s *lcp,
                         int64_t pos, int64_t *run_idx, int64_t *match_len, int8_t c);

// Compute matching statistics for pattern[0..len-1] using move + LCP.
// pattern is nt6-encoded (1-5). ms[i] = longest prefix of pattern[i..] in T.
// Returns 0 on success, -1 on error.
int rb3_move_ms_compute(const rb3_move_t *m, const struct rb3_lcp_s *lcp,
                        int64_t len, const uint8_t *pattern, int64_t *ms);

// Serialization: save move table to .mvi file; returns 0 on success, -1 on error.
int rb3_move_save(const rb3_move_t *m, const char *fn);

// Deserialization: load move table from .mvi file via memory mapping.
// Returns NULL on error. The loaded table is read-only (mmap'd).
rb3_move_t *rb3_move_load(const char *fn);

/*
 * b-move: bidirectional move structure for FMD-style bidirectional search.
 *
 * In the FMD model, both forward and backward extensions use rank queries
 * on the same BWT (exploiting symmetric construction with both strands).
 * The b-move wraps a single move structure with a persistent cumulative
 * rank table, enabling O(log r)-time rank queries at arbitrary BWT positions.
 */
typedef struct rb3_bmove_s {
	const rb3_move_t *mv;  // Move structure (not owned; caller manages lifetime)
	int64_t *cumrank;      // Cumulative rank table: cumrank[i*6+c] = rank(c, rows[i].p)
} rb3_bmove_t;

// Build b-move from existing move structure. Does not take ownership of mv.
rb3_bmove_t *rb3_bmove_init(const rb3_move_t *mv);
void rb3_bmove_destroy(rb3_bmove_t *bm);

// Rank query: compute rank(c, pos) for all characters c using b-move cumrank table.
// ok[c] = number of character c in BWT[0..pos).
void rb3_bmove_rank1a(const rb3_bmove_t *bm, int64_t pos, int64_t ok[RB3_ASIZE]);

// Dual rank query: compute rank arrays at positions k and l.
void rb3_bmove_rank2a(const rb3_bmove_t *bm, int64_t k, int64_t l, int64_t ok[RB3_ASIZE], int64_t ol[RB3_ASIZE]);

// FMD-style bidirectional extension using b-move.
// Same semantics as rb3_fmd_extend(): is_back=1 for backward, is_back=0 for forward.
void rb3_bmove_extend(const rb3_bmove_t *bm, const rb3_sai_t *ik, rb3_sai_t ok[6], int is_back);

// SMEM finding using b-move (original algorithm, same as rb3_fmd_smem).
int64_t rb3_bmove_smem(void *km, const rb3_bmove_t *bm, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);

// SMEM finding using b-move (Travis Gagie algorithm, same as rb3_fmd_smem_TG).
int64_t rb3_bmove_smem_TG(void *km, const rb3_bmove_t *bm, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);

#ifdef __cplusplus
}
#endif

#endif
