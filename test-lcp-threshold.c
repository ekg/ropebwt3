#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "lcp.h"

/*
 * Naive BWT construction: O(n^2 log n) via sorting cyclic rotations.
 * Input: nt6-encoded text with sentinel (0) at end.
 * Returns BWT as a malloc'd array.
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
 * Build FMI from a text string (nt6-encoded, not including sentinel).
 * Appends sentinel, computes BWT, builds FM-index.
 */
static int build_fmi_from_text(rb3_fmi_t *fmi, const int64_t *text, int64_t len)
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
 * Brute-force: build full SA, LCP, verify thresholds.
 */
static int verify_thresholds_brute(const rb3_fmi_t *f, const rb3_lcp_t *lcp)
{
	int64_t n = f->acc[RB3_ASIZE], i;
	int64_t *sa, *full_lcp;
	uint8_t *text;
	int64_t ok[RB3_ASIZE];
	int c, errors = 0;

	sa = (int64_t *)malloc(n * sizeof(int64_t));
	text = (uint8_t *)malloc(n);
	full_lcp = (int64_t *)calloc(n, sizeof(int64_t));

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
		int cc;
		for (cc = 0; cc < RB3_ASIZE; ++cc)
			if (i < f->acc[cc + 1]) break;
		text[sa[i]] = (uint8_t)cc;
	}

	/* Compute full LCP naively */
	for (i = 1; i < n; ++i) {
		int64_t a = sa[i - 1], b = sa[i], l = 0;
		while (a + l < n && b + l < n && text[a + l] == text[b + l])
			++l;
		full_lcp[i] = l;
	}

	/* Verify thresholds: th[i] = min(LCP[run_starts[i]], LCP[run_starts[i+1]]) */
	for (i = 0; i < lcp->n_runs; ++i) {
		int64_t left = full_lcp[lcp->run_starts[i]];
		int64_t right = (i + 1 < lcp->n_runs) ? full_lcp[lcp->run_starts[i + 1]] : 0;
		int64_t expected = left < right ? left : right;
		if (lcp->thresholds[i] != expected) {
			fprintf(stderr, "  FAIL: run %ld th=%ld expected=min(%ld,%ld)=%ld\n",
					(long)i, (long)lcp->thresholds[i], (long)left, (long)right, (long)expected);
			++errors;
		}
	}

	free(sa);
	free(text);
	free(full_lcp);
	return errors;
}

static int test_text(const char *name, const int64_t *text, int64_t len)
{
	rb3_fmi_t fmi = {0};
	rb3_lcp_t *lcp;
	int ret = 0;

	fprintf(stderr, "Test: %s (len=%ld)\n", name, (long)len);

	if (build_fmi_from_text(&fmi, text, len) < 0) {
		fprintf(stderr, "  FAIL: could not build FMI\n");
		return 1;
	}

	lcp = rb3_lcp_build(&fmi);
	if (lcp == 0) {
		fprintf(stderr, "  FAIL: rb3_lcp_build returned NULL\n");
		rb3_fmi_free(&fmi);
		return 1;
	}

	rb3_lcp_build_thresholds(lcp);

	if (lcp->thresholds == 0) {
		fprintf(stderr, "  FAIL: thresholds not allocated\n");
		ret = 1;
		goto done;
	}

	/* Brute-force verification */
	{
		int errors = verify_thresholds_brute(&fmi, lcp);
		if (errors > 0) {
			fprintf(stderr, "  FAIL: brute-force verification found %d errors\n", errors);
			ret = 1;
			goto done;
		}
	}

	/* Test rb3_lcp_threshold query matches direct access */
	{
		int64_t i;
		for (i = 0; i < lcp->n_runs; ++i) {
			int64_t th = rb3_lcp_threshold(lcp, i);
			if (th != lcp->thresholds[i]) {
				fprintf(stderr, "  FAIL: rb3_lcp_threshold(%ld)=%ld, direct=%ld\n",
						(long)i, (long)th, (long)lcp->thresholds[i]);
				ret = 1;
			}
		}
		/* Out-of-bounds queries should return 0 */
		if (rb3_lcp_threshold(lcp, -1) != 0 || rb3_lcp_threshold(lcp, lcp->n_runs) != 0) {
			fprintf(stderr, "  FAIL: out-of-bounds threshold query did not return 0\n");
			ret = 1;
		}
		if (rb3_lcp_threshold(0, 0) != 0) {
			fprintf(stderr, "  FAIL: NULL lcp threshold query did not return 0\n");
			ret = 1;
		}
	}

	/* Print summary */
	{
		int64_t i;
		fprintf(stderr, "  %ld runs, thresholds:", (long)lcp->n_runs);
		for (i = 0; i < lcp->n_runs; ++i)
			fprintf(stderr, " %ld", (long)lcp->thresholds[i]);
		fprintf(stderr, "\n");
	}

	if (ret == 0)
		fprintf(stderr, "  PASS\n");

done:
	rb3_lcp_destroy(lcp);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Test: single character A */
static int test_single_char(void)
{
	int64_t text[] = {1}; /* A */
	return test_text("A", text, 1);
}

/* Test: AAAA - all same characters */
static int test_aaaa(void)
{
	int64_t text[] = {1, 1, 1, 1}; /* AAAA */
	return test_text("AAAA", text, 4);
}

/* Test: ACAACG */
static int test_acaacg(void)
{
	int64_t text[] = {1, 2, 1, 1, 2, 3}; /* ACAACG */
	return test_text("ACAACG", text, 6);
}

/* Test: AACACA */
static int test_aacaca(void)
{
	int64_t text[] = {1, 1, 2, 1, 2, 1}; /* AACACA */
	return test_text("AACACA", text, 6);
}

/* Test: longer repetitive: ACGTACGTACGT */
static int test_repetitive(void)
{
	int64_t text[] = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4}; /* ACGTACGTACGT */
	return test_text("ACGTACGTACGT", text, 12);
}

/* Test: ACGTTTAAACCCCGGGG - varied runs */
static int test_varied(void)
{
	int64_t text[] = {1, 2, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
	return test_text("ACGTTTAAACCCCGGGG", text, 17);
}

/* Test: highly repetitive - 100 copies of AC */
static int test_highly_repetitive(void)
{
	int64_t text[200];
	int64_t i;
	for (i = 0; i < 200; i += 2) {
		text[i] = 1;     /* A */
		text[i + 1] = 2; /* C */
	}
	return test_text("(AC)x100", text, 200);
}

/* Test: build thresholds on NULL / empty */
static int test_edge_cases(void)
{
	rb3_lcp_t empty = {0};
	int ret = 0;

	fprintf(stderr, "Test: edge cases\n");

	/* NULL should not crash */
	rb3_lcp_build_thresholds(0);

	/* Empty struct should not crash */
	rb3_lcp_build_thresholds(&empty);

	/* Query on NULL */
	if (rb3_lcp_threshold(0, 0) != 0) {
		fprintf(stderr, "  FAIL: threshold on NULL didn't return 0\n");
		ret = 1;
	}

	if (ret == 0)
		fprintf(stderr, "  PASS\n");
	return ret;
}

int main(void)
{
	int ret = 0;

	rb3_verbose = 3;
	fprintf(stderr, "=== LCP Threshold Tests ===\n\n");

	ret |= test_edge_cases();
	ret |= test_single_char();
	ret |= test_aaaa();
	ret |= test_acaacg();
	ret |= test_aacaca();
	ret |= test_repetitive();
	ret |= test_varied();
	ret |= test_highly_repetitive();

	fprintf(stderr, "\n");
	if (ret == 0)
		fprintf(stderr, "=== ALL TESTS PASSED ===\n");
	else
		fprintf(stderr, "=== SOME TESTS FAILED ===\n");

	return ret;
}
