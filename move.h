#ifndef RB3_MOVE_H
#define RB3_MOVE_H

#include <stdint.h>
#include "fm-index.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Legacy 48-byte row struct, kept for v1 .mvi backward compat loading */
typedef struct rb3_move_row_s {
	int64_t p;
	int64_t pi;
	int64_t xi;
	int64_t len;
	int16_t dist[RB3_ASIZE];
	int8_t c;
	uint8_t pad_[3];
} rb3_move_row_t;

/*
 * Move table: separate arrays for cache-friendly access.
 *
 * Hot path (rank queries via bmove): binary search on p[], then scan c[]+len[]
 * to reconstruct sampled cumrank. Working set: p[] (8B/row) + c[] (1B/row) +
 * len[] (4B/row) + sampled cumrank (~0.75B/row) ≈ 14B/row vs old 96B/row.
 *
 * Cold path (LF-mapping, reposition): pi[], xi[], dist[].
 */
typedef struct rb3_move_s {
	int64_t n_runs;          /* Number of runs (r) */
	int64_t bwt_len;         /* Total BWT length (n) */
	int64_t acc[7];          /* Cumulative character counts C[] */
	int32_t d;               /* Run-splitting depth (0 = no splitting) */

	/* Separate arrays (always populated) */
	int64_t *p;              /* p[i]: starting BWT offset of run i */
	int64_t *pi;             /* pi[i]: LF[p[i]] */
	uint32_t *xi;            /* xi[i]: destination run index */
	int64_t *len;            /* len[i]: run length */
	int8_t *c;               /* c[i]: BWT character (0-5) */
	int16_t *dist;           /* dist[i*6+j]: reposition distance for char j from run i */

	int _fd;                 /* file descriptor for mmap (internal) */
	void *_mmap_base;        /* mmap base address (internal, NULL if heap-allocated) */
	size_t _mmap_len;        /* mmap region length (internal) */
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
int64_t *rb3_move_lcp_thresholds(const rb3_move_t *m, const struct rb3_lcp_s *lcp);

// Map each move row to its LCP run index.
int64_t *rb3_move_lcp_run_map(const rb3_move_t *m, const struct rb3_lcp_s *lcp);

// One backward step of matching statistics with move + LCP.
int64_t rb3_move_ms_step(const rb3_move_t *m, const int64_t *run_map,
                         const struct rb3_lcp_s *lcp,
                         int64_t pos, int64_t *run_idx, int64_t *match_len, int8_t c);

// Compute matching statistics for pattern[0..len-1] using move + LCP.
int rb3_move_ms_compute(const rb3_move_t *m, const struct rb3_lcp_s *lcp,
                        int64_t len, const uint8_t *pattern, int64_t *ms);

// Serialization: save move table to .mvi file; returns 0 on success, -1 on error.
int rb3_move_save(const rb3_move_t *m, const char *fn);

// Deserialization: load move table from .mvi file.
// Supports both v1 (48-byte rows) and v2 (compact) formats.
// Returns NULL on error.
rb3_move_t *rb3_move_load(const char *fn);

/*
 * b-move: bidirectional move structure for FMD-style bidirectional search.
 *
 * Uses a SAMPLED cumulative rank table (every CUMRANK_SAMPLE rows) to reduce
 * memory from r*48 bytes to r*48/K bytes. Rank queries reconstruct intermediate
 * values by scanning K compact rows — sequential, cache-friendly access.
 */
#define RB3_CUMRANK_SAMPLE 64

typedef struct rb3_bmove_s {
	const rb3_move_t *mv;    /* Move structure (not owned; caller manages lifetime) */
	int64_t *cumrank;        /* Sampled cumrank: cumrank[(i/K)*6+c] = rank(c, p[i*K]) */
	int64_t n_samples;       /* Number of samples: n_runs/K + 1 */
} rb3_bmove_t;

// Build b-move from existing move structure. Does not take ownership of mv.
rb3_bmove_t *rb3_bmove_init(const rb3_move_t *mv);
void rb3_bmove_destroy(rb3_bmove_t *bm);

// Rank query: compute rank(c, pos) for all characters c.
// ok[c] = number of character c in BWT[0..pos).
void rb3_bmove_rank1a(const rb3_bmove_t *bm, int64_t pos, int64_t ok[RB3_ASIZE]);

// Dual rank query: compute rank arrays at positions k and l.
void rb3_bmove_rank2a(const rb3_bmove_t *bm, int64_t k, int64_t l, int64_t ok[RB3_ASIZE], int64_t ol[RB3_ASIZE]);

// FMD-style bidirectional extension using b-move.
void rb3_bmove_extend(const rb3_bmove_t *bm, const rb3_sai_t *ik, rb3_sai_t ok[6], int is_back);

// SMEM finding using b-move (original algorithm, same as rb3_fmd_smem).
int64_t rb3_bmove_smem(void *km, const rb3_bmove_t *bm, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);

// SMEM finding using b-move (Travis Gagie algorithm, same as rb3_fmd_smem_TG).
int64_t rb3_bmove_smem_TG(void *km, const rb3_bmove_t *bm, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len);

#ifdef __cplusplus
}
#endif

#endif
