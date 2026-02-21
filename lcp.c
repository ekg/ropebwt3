#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "lcp.h"
#include "ketopt.h"

/**********************
 * Helper functions   *
 **********************/

/* Get the first-column (F) character at BWT position i */
static inline int lcp_get_f_char(const rb3_fmi_t *f, int64_t i)
{
	int c;
	for (c = 0; c < RB3_ASIZE; ++c)
		if (i < f->acc[c + 1])
			return c;
	return -1;
}

/*
 * Select: find position of the r-th occurrence (0-indexed) of character c in
 * the BWT.  Implemented via binary search on rank queries.
 */
static int64_t lcp_bwt_select(const rb3_fmi_t *f, int c, int64_t r)
{
	int64_t lo = 0, hi = f->acc[RB3_ASIZE] - 1, ok[RB3_ASIZE];
	while (lo < hi) {
		int64_t mid = lo + (hi - lo) / 2;
		rb3_fmi_rank1a(f, mid + 1, ok);
		if (ok[c] >= r + 1) hi = mid;
		else lo = mid + 1;
	}
	return lo;
}

/*
 * Compute Psi(i) = BWT position of suffix SA[i]+1.
 * Psi(i) = select_{F[i]}(i - C[F[i]]) in the BWT.
 */
static inline int64_t lcp_compute_psi(const rb3_fmi_t *f, int64_t i)
{
	int c = lcp_get_f_char(f, i);
	return lcp_bwt_select(f, c, i - f->acc[c]);
}

/*
 * Compute LCP at a run boundary at position pos (BWT[pos] != BWT[pos-1]).
 * Walk forward in both suffixes using Psi, comparing first-column characters.
 */
static int64_t lcp_at_boundary(const rb3_fmi_t *f, int64_t pos)
{
	int64_t p1 = pos - 1, p2 = pos, lcp = 0;
	while (1) {
		int c1 = lcp_get_f_char(f, p1);
		int c2 = lcp_get_f_char(f, p2);
		if (c1 != c2) break;
		if (c1 == 0) break; /* both sentinels: different sequences */
		++lcp;
		p1 = lcp_compute_psi(f, p1);
		p2 = lcp_compute_psi(f, p2);
	}
	return lcp;
}

/*
 * Compute LCP between consecutive SA entries at any position (not just
 * run boundaries). Same Psi-walk algorithm as lcp_at_boundary.
 */
static int64_t lcp_at_position(const rb3_fmi_t *f, int64_t pos)
{
	int64_t p1 = pos - 1, p2 = pos, lcp = 0;
	while (1) {
		int c1 = lcp_get_f_char(f, p1);
		int c2 = lcp_get_f_char(f, p2);
		if (c1 != c2) break;
		if (c1 == 0) break;
		++lcp;
		p1 = lcp_compute_psi(f, p1);
		p2 = lcp_compute_psi(f, p2);
	}
	return lcp;
}

/**********************
 * Public interface   *
 **********************/

rb3_lcp_t *rb3_lcp_build(const rb3_fmi_t *f)
{
	rb3_lcp_t *lcp;
	int64_t i, pos, n_runs;

	if (f->e == 0 && f->r == 0) {
		fprintf(stderr, "[E::%s] no BWT loaded\n", __func__);
		return 0;
	}

	/* Phase 1: count BWT runs */
	n_runs = rb3_fmi_get_r(f);
	if (n_runs <= 0) return 0;

	lcp = RB3_CALLOC(rb3_lcp_t, 1);
	lcp->n_runs = n_runs;
	lcp->fmi = f;
	lcp->run_starts = RB3_MALLOC(int64_t, n_runs);
	lcp->lcp_samples = RB3_MALLOC(int64_t, n_runs);

	/* Phase 2: record run boundary positions by scanning the BWT */
	if (f->e) { /* FMD backend: use rld_dec for efficient scanning */
		rlditr_t itr;
		int c;
		int64_t l;
		pos = 0; i = 0;
		rld_itr_init(f->e, &itr, 0);
		while ((l = rld_dec(f->e, &itr, &c, 0)) > 0) {
			lcp->run_starts[i++] = pos;
			pos += l;
		}
	} else { /* FMR backend: scan position by position */
		int64_t n = f->acc[RB3_ASIZE], ok[RB3_ASIZE];
		int last_c = -1, c;
		i = 0;
		for (pos = 0; pos < n; ++pos) {
			c = rb3_fmi_rank1a(f, pos, ok);
			if (c != last_c) {
				lcp->run_starts[i++] = pos;
				last_c = c;
			}
		}
	}
	assert(i == n_runs);

	/* Phase 3: compute LCP at each run boundary */
	lcp->lcp_samples[0] = 0; /* LCP[0] = 0 by convention */
	for (i = 1; i < n_runs; ++i)
		lcp->lcp_samples[i] = lcp_at_boundary(f, lcp->run_starts[i]);

	/* Phase 4: compute MONI threshold positions (tau) and within-run
	 * minimum LCP (within_min) for each run.
	 *
	 * tau[r] partitions run r into left and right zones:
	 *   - Left zone [run_starts[r], tau[r]): all within-run LCPs from
	 *     run_starts[r]+1..tau[r]-1 are >= lcp_samples[r], so the LCP to the
	 *     left boundary equals exactly lcp_samples[r].
	 *   - Right zone [tau[r], run_end): all within-run LCPs from
	 *     tau[r]+1..run_end-1 are >= lcp_samples[r+1], so the LCP to the
	 *     right boundary equals exactly lcp_samples[r+1].
	 *
	 * within_min[r] = min of all within-run LCPs (fallback for when the
	 * reposition direction doesn't match the zone). */
	{
		int64_t n = f->acc[RB3_ASIZE];
		lcp->tau = RB3_MALLOC(int64_t, n_runs);
		lcp->within_min = RB3_MALLOC(int64_t, n_runs);
		for (i = 0; i < n_runs; ++i) {
			int64_t s = lcp->run_starts[i];
			int64_t e = (i + 1 < n_runs) ? lcp->run_starts[i + 1] : n;
			int64_t right_lcp = (i + 1 < n_runs) ? lcp->lcp_samples[i + 1] : 0;
			if (e - s <= 1) {
				lcp->tau[i] = s;
				lcp->within_min[i] = INT64_MAX;
			} else {
				int64_t j, running_min, overall_min;
				/* Compute within_min and tau in a single right-to-left scan */
				lcp->tau[i] = e - 1;
				running_min = INT64_MAX;
				overall_min = INT64_MAX;
				for (j = e - 1; j > s; --j) {
					int64_t val = lcp_at_position(f, j);
					if (val < overall_min) overall_min = val;
					if (val < running_min) running_min = val;
					if (running_min >= right_lcp)
						lcp->tau[i] = j - 1;
					/* Once running_min drops below right_lcp, no more
					 * positions can join the right zone, but we continue
					 * scanning to compute overall_min */
				}
				lcp->within_min[i] = overall_min;
			}
		}
	}

	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] computed LCP at %ld run boundaries\n", __func__, (long)n_runs);
	return lcp;
}

void rb3_lcp_build_thresholds(rb3_lcp_t *lcp)
{
	int64_t i;
	if (lcp == 0 || lcp->n_runs == 0) return;
	if (lcp->thresholds) free(lcp->thresholds);
	lcp->thresholds = RB3_MALLOC(int64_t, lcp->n_runs);
	for (i = 0; i < lcp->n_runs; ++i) {
		int64_t left = lcp->lcp_samples[i];
		int64_t right = (i + 1 < lcp->n_runs) ? lcp->lcp_samples[i + 1] : 0;
		lcp->thresholds[i] = left < right ? left : right;
	}
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] computed thresholds for %ld runs\n", __func__, (long)lcp->n_runs);
}

int64_t rb3_lcp_threshold(const rb3_lcp_t *lcp, int64_t run_idx)
{
	if (lcp == 0 || lcp->thresholds == 0 || run_idx < 0 || run_idx >= lcp->n_runs) return 0;
	return lcp->thresholds[run_idx];
}

int64_t rb3_lcp_query(const rb3_lcp_t *lcp, int64_t bwt_pos)
{
	int64_t lo = 0, hi = lcp->n_runs - 1;
	/* binary search for rightmost run_starts[lo] <= bwt_pos */
	while (lo < hi) {
		int64_t mid = lo + (hi - lo + 1) / 2;
		if (lcp->run_starts[mid] <= bwt_pos) lo = mid;
		else hi = mid - 1;
	}
	return lcp->lcp_samples[lo];
}

void rb3_lcp_destroy(rb3_lcp_t *lcp)
{
	if (lcp == 0) return;
	free(lcp->run_starts);
	free(lcp->lcp_samples);
	free(lcp->thresholds);
	free(lcp->tau);
	free(lcp->within_min);
	free(lcp);
}

int64_t rb3_lcp_at_position(const rb3_lcp_t *lcp, int64_t pos)
{
	if (pos <= 0 || lcp->fmi == 0) return 0;
	return lcp_at_position(lcp->fmi, pos);
}

/***************************
 * Brute-force verification *
 ***************************/

/*
 * Build full SA and LCP arrays from the FM-index using LF walks + naive
 * character comparison.  Only for testing with small inputs.
 * Returns 0 on success, number of mismatches otherwise.
 */
static int lcp_verify(const rb3_fmi_t *f, const rb3_lcp_t *lcp)
{
	int64_t n = f->acc[RB3_ASIZE], i;
	int64_t *sa, *full_lcp;
	uint8_t *text;
	int64_t ok[RB3_ASIZE];
	int c, errors = 0;

	if (f->acc[1] != 1) {
		fprintf(stderr, "[W::%s] verification only supported for single-sequence BWTs\n", __func__);
		return -1;
	}
	if (n > 100000) {
		fprintf(stderr, "[W::%s] input too large for brute-force verification (n=%ld)\n", __func__, (long)n);
		return -1;
	}

	sa = RB3_MALLOC(int64_t, n);
	text = RB3_MALLOC(uint8_t, n);
	full_lcp = RB3_CALLOC(int64_t, n);

	/* Build SA from LF walks starting at the sentinel */
	{
		int64_t k = 0;
		for (i = n - 1; i >= 0; --i) {
			sa[k] = i;
			c = rb3_fmi_rank1a(f, k, ok);
			k = f->acc[c] + ok[c];
		}
	}

	/* Reconstruct text: T[SA[i]] = F[i] */
	for (i = 0; i < n; ++i)
		text[sa[i]] = (uint8_t)lcp_get_f_char(f, i);

	/* Compute full LCP naively */
	for (i = 1; i < n; ++i) {
		int64_t a = sa[i - 1], b = sa[i], l = 0;
		while (a + l < n && b + l < n && text[a + l] == text[b + l])
			++l;
		full_lcp[i] = l;
	}

	/* Compare with sampled values */
	for (i = 0; i < lcp->n_runs; ++i) {
		int64_t pos = lcp->run_starts[i];
		if (lcp->lcp_samples[i] != full_lcp[pos]) {
			fprintf(stderr, "[E::%s] LCP mismatch at run %ld (bwt_pos=%ld): computed=%ld expected=%ld\n",
					__func__, (long)i, (long)pos, (long)lcp->lcp_samples[i], (long)full_lcp[pos]);
			++errors;
		}
	}
	if (errors == 0)
		fprintf(stderr, "[M::%s] LCP verification passed: %ld run-boundary values correct\n", __func__, (long)lcp->n_runs);

	/* Verify thresholds if computed */
	if (lcp->thresholds) {
		int th_errors = 0;
		for (i = 0; i < lcp->n_runs; ++i) {
			int64_t left = full_lcp[lcp->run_starts[i]];
			int64_t right = (i + 1 < lcp->n_runs) ? full_lcp[lcp->run_starts[i + 1]] : 0;
			int64_t expected = left < right ? left : right;
			if (lcp->thresholds[i] != expected) {
				fprintf(stderr, "[E::%s] threshold mismatch at run %ld: computed=%ld expected=min(%ld,%ld)=%ld\n",
						__func__, (long)i, (long)lcp->thresholds[i], (long)left, (long)right, (long)expected);
				++th_errors;
			}
		}
		if (th_errors == 0)
			fprintf(stderr, "[M::%s] threshold verification passed: %ld values correct\n", __func__, (long)lcp->n_runs);
		errors += th_errors;
	}

	free(sa);
	free(text);
	free(full_lcp);
	return errors;
}

/*****************************
 * Matching statistics / PML *
 *****************************/

/* Find the run index containing BWT position pos (binary search on run_starts) */
static inline int64_t lcp_find_run(const rb3_lcp_t *lcp, int64_t pos)
{
	int64_t lo = 0, hi = lcp->n_runs - 1;
	while (lo < hi) {
		int64_t mid = lo + (hi - lo + 1) / 2;
		if (lcp->run_starts[mid] <= pos) lo = mid;
		else hi = mid - 1;
	}
	return lo;
}

/*
 * Compute matching statistics. MS[i] = length of the longest substring
 * starting at pattern[i] that occurs in the indexed text.
 *
 * Algorithm (MONI-style):
 *   Process pattern left-to-right, maintaining a BWT interval [k,l) and
 *   current match length d. For position i, try to extend the current match
 *   by prepending pattern[i] (backward search). If the interval remains
 *   non-empty, d increases. If it becomes empty, use the LCP threshold at
 *   the current run boundary to determine the new (shorter) match length,
 *   then reposition to continue.
 *
 * The pattern must be in nt6 encoding (0-5).
 */
void rb3_ms_compute(const rb3_fmi_t *f, const rb3_lcp_t *lcp, const uint8_t *pattern, int64_t len, int64_t *ms)
{
	int64_t i, k, l, d;
	int64_t ok[RB3_ASIZE], ol[RB3_ASIZE];

	if (len <= 0) return;

	/* Initialize: start with the full BWT range */
	k = 0;
	l = f->acc[RB3_ASIZE];
	d = 0;

	for (i = len - 1; i >= 0; --i) {
		int c = pattern[i];
		int64_t nk, nl;

		/* Try to extend the match by one character */
		rb3_fmi_rank2a(f, k, l, ok, ol);
		nk = f->acc[c] + ok[c];
		nl = f->acc[c] + ol[c];

		if (nk < nl) {
			/* Extension succeeded */
			k = nk;
			l = nl;
			++d;
		} else {
			/* Extension failed: shrink d using exact LCP values,
			 * widen interval, and retry.
			 *
			 * We use max(LCP[k], LCP[l]) as the shrink target: this
			 * widens the interval on at least one side (the side
			 * with the higher LCP), preserving the longest possible
			 * match.  Using min would over-shrink and miss valid
			 * intermediate matches.
			 *
			 * LCP[k] and LCP[l] are computed exactly via Psi-walk
			 * in every iteration, since widening + F-column
			 * clamping can leave k/l within a run where the
			 * run-boundary approximation would be wrong.
			 *
			 * After widening, we clamp [k, l) to the F-column
			 * character range to prevent the interval from spanning
			 * across BWT runs that straddle a character boundary. */
			while (d > 0) {
				int64_t lcp_k, lcp_l, th;
				int64_t run_idx, lo_run, hi_run;
				int fc;

				lcp_k = (k > 0) ? lcp_at_position(f, k) : 0;
				lcp_l = (l > 0 && l < f->acc[RB3_ASIZE]) ? lcp_at_position(f, l) : 0;
				th = lcp_k > lcp_l ? lcp_k : lcp_l;

				if (th < d)
					d = th;
				else
					--d; /* safety: always make progress */

				/* Find F-column character for clamping */
				for (fc = 0; fc < RB3_ASIZE; ++fc)
					if (k < f->acc[fc + 1]) break;

				/* Widen interval to all suffixes sharing a d-length prefix */
				run_idx = lcp_find_run(lcp, k);
				lo_run = run_idx;
				hi_run = lcp_find_run(lcp, l > 0 ? l - 1 : 0);
				while (lo_run > 0 && lcp->lcp_samples[lo_run] >= d)
					--lo_run;
				while (hi_run + 1 < lcp->n_runs && lcp->lcp_samples[hi_run + 1] >= d)
					++hi_run;
				k = lcp->run_starts[lo_run];
				if (hi_run + 1 < lcp->n_runs)
					l = lcp->run_starts[hi_run + 1];
				else
					l = f->acc[RB3_ASIZE];

				/* Clamp to F-column character range: prevents the
				 * interval from spanning a BWT run that straddles
				 * two different first-column characters. */
				if (d > 0) {
					if (k < f->acc[fc]) k = f->acc[fc];
					if (l > f->acc[fc + 1]) l = f->acc[fc + 1];
				}

				if (d == 0) break;

				/* Try extending again with the same character */
				rb3_fmi_rank2a(f, k, l, ok, ol);
				nk = f->acc[c] + ok[c];
				nl = f->acc[c] + ol[c];
				if (nk < nl) {
					k = nk; l = nl; ++d;
					break;
				}
			}
			if (d == 0) {
				k = f->acc[c]; l = f->acc[c + 1];
				if (k < l) d = 1;
			}
		}
		ms[i] = d;
	}
}

/*
 * Compute pseudo-matching lengths (PML). Similar to MS but simpler:
 * when a match cannot be extended, record the current length as PML,
 * then reposition to continue from a shorter match.
 *
 * PML[i] <= MS[i] for all i. PML is faster because we don't need to
 * find the exact longest match -- we just take whatever the threshold gives us.
 *
 * The pattern must be in nt6 encoding (0-5).
 */
void rb3_pml_compute(const rb3_fmi_t *f, const rb3_lcp_t *lcp, const uint8_t *pattern, int64_t len, int64_t *pml)
{
	int64_t i, k, l, d;
	int64_t ok[RB3_ASIZE], ol[RB3_ASIZE];

	if (len <= 0) return;

	k = 0;
	l = f->acc[RB3_ASIZE];
	d = 0;

	for (i = len - 1; i >= 0; --i) {
		int c = pattern[i];
		int64_t nk, nl;

		rb3_fmi_rank2a(f, k, l, ok, ol);
		nk = f->acc[c] + ok[c];
		nl = f->acc[c] + ol[c];

		if (nk < nl) {
			k = nk;
			l = nl;
			++d;
		} else {
			/* Record current match length, then shrink using threshold */
			int64_t run_idx = lcp_find_run(lcp, k);
			int64_t th = rb3_lcp_threshold(lcp, run_idx);
			d = th < d ? th : d;

			if (d > 0) {
				/* Widen interval to match the shorter prefix */
				int64_t lo_run = run_idx, hi_run = run_idx;
				while (lo_run > 0 && lcp->lcp_samples[lo_run] >= d)
					--lo_run;
				while (hi_run + 1 < lcp->n_runs && lcp->lcp_samples[hi_run + 1] >= d)
					++hi_run;
				k = lcp->run_starts[lo_run];
				l = (hi_run + 1 < lcp->n_runs) ? lcp->run_starts[hi_run + 1] : f->acc[RB3_ASIZE];

				/* Try extending with current character */
				rb3_fmi_rank2a(f, k, l, ok, ol);
				nk = f->acc[c] + ok[c];
				nl = f->acc[c] + ol[c];
				if (nk < nl) {
					k = nk;
					l = nl;
					++d;
				}
				/* If still can't extend, keep d as is (PML approximation) */
			}
			if (d == 0) {
				k = f->acc[c];
				l = f->acc[c + 1];
				if (k < l) d = 1;
			}
		}
		pml[i] = d;
	}
}

/*******************
 * main() function *
 *******************/

int main_lcp(int argc, char *argv[])
{
	int c, verify = 0, do_thresholds = 0;
	rb3_fmi_t f;
	rb3_lcp_t *lcp;
	ketopt_t o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, "vt", 0)) >= 0) {
		if (c == 'v') verify = 1;
		else if (c == 't') do_thresholds = 1;
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: ropebwt3 lcp [options] <in.fmd>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t       compute thresholds for matching statistics\n");
		fprintf(stderr, "  -v       verify against brute-force (small inputs only)\n");
		return 1;
	}

	rb3_fmi_restore(&f, argv[o.ind], 0);
	if (f.e == 0 && f.r == 0) {
		fprintf(stderr, "[E::%s] failed to load the FM-index\n", __func__);
		return 1;
	}

	lcp = rb3_lcp_build(&f);
	if (lcp == 0) {
		fprintf(stderr, "[E::%s] failed to build LCP\n", __func__);
		rb3_fmi_free(&f);
		return 1;
	}

	if (do_thresholds)
		rb3_lcp_build_thresholds(lcp);

	/* Print LCP values (and thresholds if computed) at run boundaries */
	{
		int64_t i;
		printf("n_runs\t%ld\n", (long)lcp->n_runs);
		for (i = 0; i < lcp->n_runs; ++i) {
			printf("%ld\t%ld", (long)lcp->run_starts[i], (long)lcp->lcp_samples[i]);
			if (lcp->thresholds)
				printf("\t%ld", (long)lcp->thresholds[i]);
			printf("\n");
		}
	}

	if (verify) {
		int ret = lcp_verify(&f, lcp);
		if (ret > 0) {
			fprintf(stderr, "[E::%s] verification FAILED with %d errors\n", __func__, ret);
			rb3_lcp_destroy(lcp);
			rb3_fmi_free(&f);
			return 1;
		}
	}

	rb3_lcp_destroy(lcp);
	rb3_fmi_free(&f);
	return 0;
}
