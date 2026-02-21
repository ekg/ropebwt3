#ifndef RB3_LCP_H
#define RB3_LCP_H

#include <stdint.h>
#include "fm-index.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rb3_lcp_s {
	int64_t n_runs;       /* number of BWT runs */
	int64_t *lcp_samples; /* LCP values at run boundaries, size n_runs */
	int64_t *run_starts;  /* BWT positions of run boundaries, size n_runs */
	int64_t *thresholds;  /* threshold values at run boundaries, size n_runs */
	int64_t *tau;         /* MONI threshold positions per run, size n_runs */
	int64_t *within_min;  /* minimum within-run LCP per run, size n_runs */
	const rb3_fmi_t *fmi; /* back-pointer to FM-index (for on-the-fly LCP) */
} rb3_lcp_t;

rb3_lcp_t *rb3_lcp_build(const rb3_fmi_t *f);
int64_t rb3_lcp_query(const rb3_lcp_t *lcp, int64_t bwt_pos);
void rb3_lcp_build_thresholds(rb3_lcp_t *lcp);
int64_t rb3_lcp_threshold(const rb3_lcp_t *lcp, int64_t run_idx);
void rb3_lcp_destroy(rb3_lcp_t *lcp);

/* Compute LCP between consecutive SA entries at arbitrary BWT position pos.
 * Returns LCP[pos] = length of longest common prefix of SA[pos-1] and SA[pos].
 * Uses the stored FM-index pointer for Psi-walk computation. */
int64_t rb3_lcp_at_position(const rb3_lcp_t *lcp, int64_t pos);

/* Matching statistics: MS[i] = length of longest match starting at pattern[i] in the indexed text */
void rb3_ms_compute(const rb3_fmi_t *f, const rb3_lcp_t *lcp, const uint8_t *pattern, int64_t len, int64_t *ms);

/* Pseudo-matching lengths: faster but less precise variant of matching statistics */
void rb3_pml_compute(const rb3_fmi_t *f, const rb3_lcp_t *lcp, const uint8_t *pattern, int64_t len, int64_t *pml);

int main_lcp(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
