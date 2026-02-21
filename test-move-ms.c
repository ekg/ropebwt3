#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "move.h"
#include "lcp.h"

/*
 * Naive BWT construction: O(n^2 log n) via sorting cyclic rotations.
 * Input: nt6-encoded text with sentinel (0) at end.
 */
static int64_t *g_text;
static int64_t g_n;

static int cmp_rotations(const void *a_, const void *b_)
{
	int64_t a = *(const int64_t *)a_, b = *(const int64_t *)b_;
	int64_t i;
	for (i = 0; i < g_n; ++i) {
		int ca = g_text[(a + i) % g_n], cb = g_text[(b + i) % g_n];
		if (ca != cb) return ca < cb ? -1 : 1;
	}
	return 0;
}

static uint8_t *naive_bwt(const int64_t *text, int64_t n)
{
	int64_t *sa, i;
	uint8_t *bwt;

	sa = (int64_t *)malloc(n * sizeof(int64_t));
	bwt = (uint8_t *)malloc(n);
	for (i = 0; i < n; ++i) sa[i] = i;

	g_text = (int64_t *)text;
	g_n = n;
	qsort(sa, n, sizeof(int64_t), cmp_rotations);

	for (i = 0; i < n; ++i)
		bwt[i] = (uint8_t)text[(sa[i] + n - 1) % n];

	free(sa);
	return bwt;
}

/*
 * Build an FM-index from nt6-encoded text (no sentinel).
 * Appends sentinel, computes BWT, builds FMD.
 */
static int build_fmi(rb3_fmi_t *fmi, const int64_t *text, int64_t len)
{
	int64_t *full_text, n = len + 1;
	uint8_t *bwt;
	rld_t *e;

	full_text = (int64_t *)malloc(n * sizeof(int64_t));
	memcpy(full_text, text, len * sizeof(int64_t));
	full_text[len] = 0; /* sentinel */

	bwt = naive_bwt(full_text, n);
	free(full_text);

	e = rb3_enc_plain2rld(n, bwt, 3);
	free(bwt);
	if (e == 0) return -1;

	rb3_fmi_init(fmi, e, 0);
	return 0;
}

/*
 * Brute-force matching statistics:
 * For each position i in pattern, find the longest prefix of pattern[i..]
 * that occurs as a substring of text. Uses text reconstruction + naive search.
 */
static void brute_ms(const rb3_fmi_t *f, int64_t plen, const uint8_t *pattern, int64_t *ms)
{
	int64_t n = f->acc[RB3_ASIZE], i, t;
	int64_t *sa;
	uint8_t *text;
	int64_t ok[RB3_ASIZE];
	int c;

	sa = (int64_t *)malloc(n * sizeof(int64_t));
	text = (uint8_t *)malloc(n);

	/* Build SA via LF walks from the sentinel */
	{
		int64_t k = 0;
		for (i = n - 1; i >= 0; --i) {
			sa[k] = i;
			c = rb3_fmi_rank1a(f, k, ok);
			k = f->acc[c] + ok[c];
		}
	}

	/* Reconstruct text: T[SA[i]] = F[i] */
	for (i = 0; i < n; ++i) {
		for (c = 0; c < RB3_ASIZE; ++c)
			if (i < f->acc[c + 1]) break;
		text[sa[i]] = (uint8_t)c;
	}

	/* For each position i in pattern, find longest prefix in text */
	for (i = 0; i < plen; ++i) {
		int64_t max_len = 0;
		for (t = 0; t < n; ++t) {
			int64_t l = 0;
			while (i + l < plen && t + l < n &&
			       pattern[i + l] == text[t + l] &&
			       pattern[i + l] > 0 && text[t + l] > 0) /* don't match sentinels */
				++l;
			if (l > max_len) max_len = l;
		}
		ms[i] = max_len;
	}

	free(sa);
	free(text);
}

/*
 * Test helper: build move+LCP, compute MS, compare with brute-force.
 * Returns 0 on success, 1 on failure.
 */
static int test_ms_text_pattern(const char *name, const int64_t *text, int64_t tlen,
                                const uint8_t *pattern, int64_t plen, int split_d)
{
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_lcp_t *lcp;
	int64_t *ms_move, *ms_brute;
	int ret = 0;
	int64_t i;

	fprintf(stderr, "Test: %s (tlen=%ld, plen=%ld, d=%d)\n",
	        name, (long)tlen, (long)plen, split_d);

	if (build_fmi(&fmi, text, tlen) < 0) {
		fprintf(stderr, "  FAIL: could not build FMI\n");
		return 1;
	}

	/* Build move table */
	m = rb3_move_build(&fmi);
	if (split_d > 0) rb3_move_split(m, split_d);
	rb3_move_precompute_dist(m);

	/* Build LCP + thresholds */
	lcp = rb3_lcp_build(&fmi);
	if (lcp == 0) {
		fprintf(stderr, "  FAIL: rb3_lcp_build returned NULL\n");
		rb3_move_destroy(m);
		rb3_fmi_free(&fmi);
		return 1;
	}
	rb3_lcp_build_thresholds(lcp);

	/* Compute MS via move+LCP */
	ms_move = (int64_t *)malloc(plen * sizeof(int64_t));
	if (rb3_move_ms_compute(m, lcp, plen, pattern, ms_move) != 0) {
		fprintf(stderr, "  FAIL: rb3_move_ms_compute returned error\n");
		ret = 1;
		goto done;
	}

	/* Compute MS via brute-force */
	ms_brute = (int64_t *)malloc(plen * sizeof(int64_t));
	brute_ms(&fmi, plen, pattern, ms_brute);

	/* Compare */
	for (i = 0; i < plen; ++i) {
		if (ms_move[i] != ms_brute[i]) {
			fprintf(stderr, "  FAIL: position %ld: move=%ld, brute=%ld\n",
			        (long)i, (long)ms_move[i], (long)ms_brute[i]);
			/* Print context */
			fprintf(stderr, "  pattern[%ld..] =", (long)i);
			{
				int64_t j;
				for (j = i; j < plen && j < i + 10; ++j)
					fprintf(stderr, " %d", pattern[j]);
			}
			fprintf(stderr, "\n");
			ret = 1;
			goto done2;
		}
	}

	/* Print summary */
	{
		fprintf(stderr, "  MS:");
		for (i = 0; i < plen && i < 20; ++i)
			fprintf(stderr, " %ld", (long)ms_move[i]);
		if (plen > 20) fprintf(stderr, " ...");
		fprintf(stderr, "\n");
	}

	fprintf(stderr, "  PASS\n");

done2:
	free(ms_brute);
done:
	free(ms_move);
	rb3_lcp_destroy(lcp);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Test: exact substring match */
static int test_exact_match(void)
{
	/* Text: ACGT (nt6: 1,2,3,4) */
	int64_t text[] = {1, 2, 3, 4};
	/* Pattern: CGT (nt6: 2,3,4) — exact substring */
	uint8_t pat[] = {2, 3, 4};
	return test_ms_text_pattern("exact_match_CGT_in_ACGT", text, 4, pat, 3, 0);
}

/* Test: partial match */
static int test_partial_match(void)
{
	/* Text: ACGT */
	int64_t text[] = {1, 2, 3, 4};
	/* Pattern: ACGA — matches ACG (len 3) at start, then A (len 1) at end */
	uint8_t pat[] = {1, 2, 3, 1};
	return test_ms_text_pattern("partial_ACGA_in_ACGT", text, 4, pat, 4, 0);
}

/* Test: no match (character not in text) */
static int test_no_match(void)
{
	/* Text: AAA (nt6: 1,1,1) */
	int64_t text[] = {1, 1, 1};
	/* Pattern: CCC (nt6: 2,2,2) — C doesn't appear as a non-sentinel */
	/* Wait, C might appear in the BWT. Let me use N (5) instead */
	uint8_t pat[] = {5, 5, 5};
	return test_ms_text_pattern("no_match_NNN_in_AAA", text, 3, pat, 3, 0);
}

/* Test: single character */
static int test_single_char(void)
{
	int64_t text[] = {1, 2, 3, 4};
	uint8_t pat[] = {2}; /* C */
	return test_ms_text_pattern("single_C_in_ACGT", text, 4, pat, 1, 0);
}

/* Test: repetitive text */
static int test_repetitive(void)
{
	/* Text: ACACAC */
	int64_t text[] = {1, 2, 1, 2, 1, 2};
	/* Pattern: ACAC — should match fully */
	uint8_t pat[] = {1, 2, 1, 2};
	return test_ms_text_pattern("ACAC_in_ACACAC", text, 6, pat, 4, 0);
}

/* Test: pattern longer than matching substring */
static int test_longer_pattern(void)
{
	/* Text: AC */
	int64_t text[] = {1, 2};
	/* Pattern: ACGT — AC matches (2), G matches (1), T matches (1) */
	uint8_t pat[] = {1, 2, 3, 4};
	return test_ms_text_pattern("ACGT_in_AC", text, 2, pat, 4, 0);
}

/* Test: with run splitting (d=2) */
static int test_split_d2(void)
{
	/* Text: ACGTACGTACGT */
	int64_t text[] = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
	uint8_t pat[] = {1, 2, 3, 4, 1, 2};
	return test_ms_text_pattern("split_d2_ACGTAC_in_ACGT3", text, 12, pat, 6, 2);
}

/* Test: split vs unsplit give same results */
static int test_split_consistency(void)
{
	int64_t text[] = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
	uint8_t pat[] = {2, 3, 4, 1, 2, 3};
	rb3_fmi_t fmi = {0};
	rb3_move_t *m0, *m2;
	rb3_lcp_t *lcp;
	int64_t *ms0, *ms2;
	int64_t plen = 6, i;
	int ret = 0;

	fprintf(stderr, "Test: split_consistency\n");

	if (build_fmi(&fmi, text, 12) < 0) {
		fprintf(stderr, "  FAIL: build_fmi\n");
		return 1;
	}

	lcp = rb3_lcp_build(&fmi);
	rb3_lcp_build_thresholds(lcp);

	/* Unsplit */
	m0 = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m0);
	ms0 = (int64_t *)malloc(plen * sizeof(int64_t));
	rb3_move_ms_compute(m0, lcp, plen, pat, ms0);

	/* Split d=2 */
	m2 = rb3_move_build(&fmi);
	rb3_move_split(m2, 2);
	rb3_move_precompute_dist(m2);
	ms2 = (int64_t *)malloc(plen * sizeof(int64_t));
	rb3_move_ms_compute(m2, lcp, plen, pat, ms2);

	for (i = 0; i < plen; ++i) {
		if (ms0[i] != ms2[i]) {
			fprintf(stderr, "  FAIL: pos %ld: unsplit=%ld, split=%ld\n",
			        (long)i, (long)ms0[i], (long)ms2[i]);
			ret = 1;
			break;
		}
	}

	if (ret == 0) fprintf(stderr, "  PASS\n");

	free(ms0);
	free(ms2);
	rb3_move_destroy(m0);
	rb3_move_destroy(m2);
	rb3_lcp_destroy(lcp);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Test: all single-character patterns against a varied text */
static int test_all_single_chars(void)
{
	int64_t text[] = {1, 2, 3, 4, 1, 1, 2, 2, 3, 3, 4, 4};
	int c, ret = 0;

	for (c = 1; c <= 5; ++c) {
		uint8_t pat[1] = {(uint8_t)c};
		char name[32];
		snprintf(name, sizeof(name), "single_char_%d", c);
		ret |= test_ms_text_pattern(name, text, 12, pat, 1, 0);
	}
	return ret;
}

/* Test: exhaustive short patterns on a small text */
static int test_exhaustive_short(void)
{
	int64_t text[] = {1, 2, 1, 1, 2, 3}; /* ACAACG */
	int c1, c2, c3, ret = 0, n_tested = 0;

	fprintf(stderr, "Test: exhaustive_short (all 2-3 char patterns on ACAACG)\n");

	/* All 2-character patterns */
	for (c1 = 1; c1 <= 4; ++c1) {
		for (c2 = 1; c2 <= 4; ++c2) {
			rb3_fmi_t fmi = {0};
			rb3_move_t *m;
			rb3_lcp_t *lcp;
			int64_t ms_move[2], ms_brute[2];
			uint8_t pat[2] = {(uint8_t)c1, (uint8_t)c2};

			if (build_fmi(&fmi, text, 6) < 0) return 1;
			m = rb3_move_build(&fmi);
			rb3_move_precompute_dist(m);
			lcp = rb3_lcp_build(&fmi);
			rb3_lcp_build_thresholds(lcp);

			rb3_move_ms_compute(m, lcp, 2, pat, ms_move);
			brute_ms(&fmi, 2, pat, ms_brute);

			if (ms_move[0] != ms_brute[0] || ms_move[1] != ms_brute[1]) {
				fprintf(stderr, "  FAIL: pat=[%d,%d] move=[%ld,%ld] brute=[%ld,%ld]\n",
				        c1, c2, (long)ms_move[0], (long)ms_move[1],
				        (long)ms_brute[0], (long)ms_brute[1]);
				ret = 1;
			}
			n_tested++;

			rb3_lcp_destroy(lcp);
			rb3_move_destroy(m);
			rb3_fmi_free(&fmi);
			if (ret) return ret;
		}
	}

	/* All 3-character patterns */
	for (c1 = 1; c1 <= 4; ++c1) {
		for (c2 = 1; c2 <= 4; ++c2) {
			for (c3 = 1; c3 <= 4; ++c3) {
				rb3_fmi_t fmi = {0};
				rb3_move_t *m;
				rb3_lcp_t *lcp;
				int64_t ms_move[3], ms_brute[3];
				uint8_t pat[3] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3};
				int i;

				if (build_fmi(&fmi, text, 6) < 0) return 1;
				m = rb3_move_build(&fmi);
				rb3_move_precompute_dist(m);
				lcp = rb3_lcp_build(&fmi);
				rb3_lcp_build_thresholds(lcp);

				rb3_move_ms_compute(m, lcp, 3, pat, ms_move);
				brute_ms(&fmi, 3, pat, ms_brute);

				for (i = 0; i < 3; ++i) {
					if (ms_move[i] != ms_brute[i]) {
						fprintf(stderr, "  FAIL: pat=[%d,%d,%d] pos=%d move=%ld brute=%ld\n",
						        c1, c2, c3, i, (long)ms_move[i], (long)ms_brute[i]);
						ret = 1;
						break;
					}
				}
				n_tested++;

				rb3_lcp_destroy(lcp);
				rb3_move_destroy(m);
				rb3_fmi_free(&fmi);
				if (ret) return ret;
			}
		}
	}

	if (ret == 0)
		fprintf(stderr, "  PASS (%d patterns tested)\n", n_tested);
	return ret;
}

/* Test: longer repetitive text with longer pattern */
static int test_long_repetitive(void)
{
	int64_t text[100], i;
	uint8_t pat[30];

	/* Text: (ACGT)x25 */
	for (i = 0; i < 100; ++i)
		text[i] = (i % 4) + 1;

	/* Pattern: ACGTACGTACGTAC (14 chars) — should all be long matches */
	for (i = 0; i < 14; ++i)
		pat[i] = (i % 4) + 1;

	return test_ms_text_pattern("long_rep_ACGT25", text, 100, pat, 14, 0);
}

/* Test: long repetitive with splitting */
static int test_long_repetitive_split(void)
{
	int64_t text[100], i;
	uint8_t pat[20];

	for (i = 0; i < 100; ++i)
		text[i] = (i % 4) + 1;

	for (i = 0; i < 20; ++i)
		pat[i] = (i % 4) + 1;

	return test_ms_text_pattern("long_rep_split_d2", text, 100, pat, 20, 2);
}

/* Test: threshold precomputation */
static int test_threshold_precompute(void)
{
	int64_t text[] = {1, 2, 1, 1, 2, 3}; /* ACAACG */
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_lcp_t *lcp;
	int64_t *th;
	int64_t i;
	int ret = 0;

	fprintf(stderr, "Test: threshold_precompute\n");

	if (build_fmi(&fmi, text, 6) < 0) return 1;
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);
	lcp = rb3_lcp_build(&fmi);
	rb3_lcp_build_thresholds(lcp);

	th = rb3_move_lcp_thresholds(m, lcp);
	if (th == 0) {
		fprintf(stderr, "  FAIL: rb3_move_lcp_thresholds returned NULL\n");
		ret = 1;
		goto done;
	}

	/* Verify: each threshold should match the LCP run containing the move row */
	for (i = 0; i < m->n_runs; ++i) {
		/* Find the LCP run for this move row's position */
		int64_t j, lcp_run = 0;
		for (j = 0; j < lcp->n_runs; ++j) {
			if (j + 1 < lcp->n_runs && lcp->run_starts[j + 1] <= m->rows[i].p)
				continue;
			lcp_run = j;
			break;
		}
		if (th[i] != lcp->thresholds[lcp_run]) {
			fprintf(stderr, "  FAIL: move row %ld (p=%ld): th=%ld, expected=%ld (lcp run %ld)\n",
			        (long)i, (long)m->rows[i].p, (long)th[i],
			        (long)lcp->thresholds[lcp_run], (long)lcp_run);
			ret = 1;
			break;
		}
	}

	if (ret == 0) fprintf(stderr, "  PASS\n");

	free(th);
done:
	rb3_lcp_destroy(lcp);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Test: error handling */
static int test_error_handling(void)
{
	int64_t text[] = {1, 2, 3};
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_lcp_t *lcp;
	int64_t ms[3];
	uint8_t pat[] = {1, 2, 3};
	int ret = 0;

	fprintf(stderr, "Test: error_handling\n");

	if (build_fmi(&fmi, text, 3) < 0) return 1;
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);
	lcp = rb3_lcp_build(&fmi);
	/* Don't build thresholds — should cause error */

	if (rb3_move_ms_compute(m, lcp, 3, pat, ms) != -1) {
		fprintf(stderr, "  FAIL: expected error for missing thresholds\n");
		ret = 1;
	}

	/* NULL move table */
	rb3_lcp_build_thresholds(lcp);
	if (rb3_move_ms_compute(NULL, lcp, 3, pat, ms) != -1) {
		fprintf(stderr, "  FAIL: expected error for NULL move\n");
		ret = 1;
	}

	/* NULL LCP */
	if (rb3_move_ms_compute(m, NULL, 3, pat, ms) != -1) {
		fprintf(stderr, "  FAIL: expected error for NULL lcp\n");
		ret = 1;
	}

	/* Zero-length pattern (should succeed) */
	if (rb3_move_ms_compute(m, lcp, 0, pat, ms) != 0) {
		fprintf(stderr, "  FAIL: zero-length should return 0\n");
		ret = 1;
	}

	if (ret == 0) fprintf(stderr, "  PASS\n");

	rb3_lcp_destroy(lcp);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Test: varied text with all character types */
static int test_varied_text(void)
{
	/* ACGTTTAAACCCCGGGG */
	int64_t text[] = {1, 2, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
	uint8_t pat[] = {4, 4, 1, 1, 2, 2, 3, 3};
	return test_ms_text_pattern("varied_text", text, 17, pat, 8, 0);
}

/* Test: pattern with mismatches forcing threshold use */
static int test_threshold_use(void)
{
	/* Text: AACACA */
	int64_t text[] = {1, 1, 2, 1, 2, 1};
	/* Pattern: AACGT — first 3 chars match (AAC), then mismatch */
	uint8_t pat[] = {1, 1, 2, 3, 4};
	return test_ms_text_pattern("threshold_AACGT_in_AACACA", text, 6, pat, 5, 0);
}

/* Test: ms_step function directly */
static int test_ms_step_direct(void)
{
	int64_t text[] = {1, 2, 3, 4};
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_lcp_t *lcp;
	int64_t *run_map;
	int64_t pos, run_idx, match_len, new_pos;
	int ret = 0;

	fprintf(stderr, "Test: ms_step_direct\n");

	if (build_fmi(&fmi, text, 4) < 0) return 1;
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);
	lcp = rb3_lcp_build(&fmi);
	rb3_lcp_build_thresholds(lcp);
	run_map = rb3_move_lcp_run_map(m, lcp);

	/* Initial state */
	pos = 0; run_idx = 0; match_len = 0;

	/* Step with character that exists */
	new_pos = rb3_move_ms_step(m, run_map, lcp, pos, &run_idx, &match_len, 1);
	if (new_pos < 0) {
		fprintf(stderr, "  FAIL: ms_step returned -1 for existing char\n");
		ret = 1;
		goto done;
	}
	if (match_len != 1) {
		fprintf(stderr, "  FAIL: first step match_len=%ld, expected 1\n", (long)match_len);
		ret = 1;
		goto done;
	}

	/* Step with non-existent character (N=5 if not in text) */
	{
		int64_t p2 = pos, r2 = run_idx, ml2 = match_len;
		int64_t res = rb3_move_ms_step(m, run_map, lcp, p2, &r2, &ml2, 5);
		if (res != -1 && fmi.acc[5] == fmi.acc[6]) {
			/* N doesn't exist, should return -1 */
			fprintf(stderr, "  FAIL: ms_step should return -1 for non-existent char\n");
			ret = 1;
			goto done;
		}
	}

	if (ret == 0) fprintf(stderr, "  PASS\n");

done:
	free(run_map);
	rb3_lcp_destroy(lcp);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Debug: trace the varied_text failure step by step */
static int test_varied_text_debug(void)
{
	/* ACGTTTAAACCCCGGGG */
	int64_t text[] = {1, 2, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
	uint8_t pat[] = {4, 4, 1, 1, 2, 2, 3, 3};
	int64_t tlen = 17, plen = 8;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_lcp_t *lcp;
	int64_t *run_map, *ms_bwt, *ms_brute;
	int64_t i, pos, run_idx, match_len;

	fprintf(stderr, "\nTest: varied_text_debug\n");

	if (build_fmi(&fmi, text, tlen) < 0) return 1;
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);
	lcp = rb3_lcp_build(&fmi);
	rb3_lcp_build_thresholds(lcp);
	run_map = rb3_move_lcp_run_map(m, lcp);

	/* Dump LCP run info */
	fprintf(stderr, "  LCP runs: %ld\n", (long)lcp->n_runs);
	for (i = 0; i < lcp->n_runs; ++i)
		fprintf(stderr, "    lcp_run[%ld]: start=%ld lcp=%ld th=%ld\n",
		        (long)i, (long)lcp->run_starts[i], (long)lcp->lcp_samples[i],
		        (long)lcp->thresholds[i]);

	/* Dump move rows */
	fprintf(stderr, "  Move rows: %ld\n", (long)m->n_runs);
	for (i = 0; i < m->n_runs; ++i)
		fprintf(stderr, "    row[%ld]: c=%d p=%ld len=%ld lcp_run=%ld dist_A=%d dist_C=%d dist_G=%d dist_T=%d\n",
		        (long)i, m->rows[i].c, (long)m->rows[i].p, (long)m->rows[i].len,
		        (long)run_map[i],
		        m->rows[i].dist[1], m->rows[i].dist[2],
		        m->rows[i].dist[3], m->rows[i].dist[4]);

	/* Trace move-based MS step by step */
	pos = 0; run_idx = 0; match_len = 0;
	fprintf(stderr, "  Move trace:\n");
	for (i = plen - 1; i >= 0; --i) {
		int8_t c = pat[i];
		int64_t old_pos = pos, old_run = run_idx, old_ml = match_len;
		pos = rb3_move_ms_step(m, run_map, lcp, pos, &run_idx, &match_len, c);
		fprintf(stderr, "    i=%ld c=%d: pos %ld(run %ld,c=%d) ml=%ld -> pos %ld(run %ld) ml=%ld\n",
		        (long)i, c, (long)old_pos, (long)old_run, m->rows[old_run].c,
		        (long)old_ml, (long)pos, (long)run_idx, (long)match_len);
	}

	/* Compare with BWT-based MS */
	ms_bwt = (int64_t *)malloc(plen * sizeof(int64_t));
	rb3_ms_compute(&fmi, lcp, pat, plen, ms_bwt);
	ms_brute = (int64_t *)malloc(plen * sizeof(int64_t));
	brute_ms(&fmi, plen, pat, ms_brute);

	fprintf(stderr, "  BWT MS: ");
	for (i = 0; i < plen; ++i) fprintf(stderr, " %ld", (long)ms_bwt[i]);
	fprintf(stderr, "\n  Brute MS:");
	for (i = 0; i < plen; ++i) fprintf(stderr, " %ld", (long)ms_brute[i]);
	fprintf(stderr, "\n");

	free(ms_bwt);
	free(ms_brute);
	free(run_map);
	rb3_lcp_destroy(lcp);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return 0;
}

int main(void)
{
	int ret = 0;

	rb3_verbose = 3;
	fprintf(stderr, "=== Move+LCP Matching Statistics Tests ===\n\n");

	ret |= test_error_handling();
	ret |= test_threshold_precompute();
	ret |= test_ms_step_direct();
	ret |= test_single_char();
	ret |= test_exact_match();
	ret |= test_partial_match();
	ret |= test_no_match();
	ret |= test_repetitive();
	ret |= test_longer_pattern();
	ret |= test_varied_text();
	ret |= test_threshold_use();
	ret |= test_all_single_chars();
	ret |= test_split_d2();
	ret |= test_split_consistency();
	ret |= test_long_repetitive();
	ret |= test_long_repetitive_split();
	ret |= test_exhaustive_short();
	test_varied_text_debug();

	fprintf(stderr, "\n");
	if (ret == 0)
		fprintf(stderr, "=== ALL TESTS PASSED ===\n");
	else
		fprintf(stderr, "=== SOME TESTS FAILED ===\n");

	return ret;
}
