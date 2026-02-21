#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "lcp.h"

/*
 * Naive BWT construction for testing: sort cyclic rotations.
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
 * Brute-force matching statistics: for each position i in pattern,
 * find the longest substring starting at i that occurs in text.
 */
static void naive_ms(const uint8_t *text, int64_t tlen,
                     const uint8_t *pattern, int64_t plen,
                     int64_t *ms)
{
	int64_t i, j, k;
	for (i = 0; i < plen; ++i) {
		ms[i] = 0;
		/* Try all starting positions in text */
		for (j = 0; j < tlen; ++j) {
			int64_t match = 0;
			for (k = 0; i + k < plen && j + k < tlen; ++k) {
				if (pattern[i + k] != text[j + k]) break;
				++match;
			}
			if (match > ms[i]) ms[i] = match;
		}
	}
}

/*
 * Test MS computation against brute-force.
 */
static int test_ms(const char *name, const int64_t *text_nt6, int64_t tlen,
                   const uint8_t *pattern_nt6, int64_t plen)
{
	rb3_fmi_t fmi = {0};
	rb3_lcp_t *lcp;
	int64_t *ms_computed, *ms_expected;
	uint8_t *text_u8;
	int64_t i;
	int ret = 0;

	fprintf(stderr, "Test MS: %s (tlen=%ld, plen=%ld)\n", name, (long)tlen, (long)plen);

	if (build_fmi_from_text(&fmi, text_nt6, tlen) < 0) {
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

	/* Compute MS using our function */
	ms_computed = (int64_t *)malloc(plen * sizeof(int64_t));
	rb3_ms_compute(&fmi, lcp, pattern_nt6, plen, ms_computed);

	/* Compute MS brute-force */
	/* Convert text from nt6 int64_t to uint8_t for comparison */
	text_u8 = (uint8_t *)malloc(tlen);
	for (i = 0; i < tlen; ++i) text_u8[i] = (uint8_t)text_nt6[i];

	ms_expected = (int64_t *)malloc(plen * sizeof(int64_t));
	naive_ms(text_u8, tlen, pattern_nt6, plen, ms_expected);

	/* Compare */
	for (i = 0; i < plen; ++i) {
		if (ms_computed[i] != ms_expected[i]) {
			fprintf(stderr, "  FAIL: MS[%ld] = %ld, expected %ld\n",
					(long)i, (long)ms_computed[i], (long)ms_expected[i]);
			ret = 1;
		}
	}

	if (ret == 0) {
		fprintf(stderr, "  MS values:");
		for (i = 0; i < plen; ++i) fprintf(stderr, " %ld", (long)ms_computed[i]);
		fprintf(stderr, "\n  PASS\n");
	}

	free(ms_computed);
	free(ms_expected);
	free(text_u8);
	rb3_lcp_destroy(lcp);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test PML: verify PML[i] <= MS[i] for all i.
 */
static int test_pml(const char *name, const int64_t *text_nt6, int64_t tlen,
                    const uint8_t *pattern_nt6, int64_t plen)
{
	rb3_fmi_t fmi = {0};
	rb3_lcp_t *lcp;
	int64_t *ms, *pml;
	int64_t i;
	int ret = 0;

	fprintf(stderr, "Test PML: %s (tlen=%ld, plen=%ld)\n", name, (long)tlen, (long)plen);

	if (build_fmi_from_text(&fmi, text_nt6, tlen) < 0) {
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

	ms = (int64_t *)malloc(plen * sizeof(int64_t));
	pml = (int64_t *)malloc(plen * sizeof(int64_t));

	rb3_ms_compute(&fmi, lcp, pattern_nt6, plen, ms);
	rb3_pml_compute(&fmi, lcp, pattern_nt6, plen, pml);

	for (i = 0; i < plen; ++i) {
		if (pml[i] > ms[i]) {
			fprintf(stderr, "  FAIL: PML[%ld]=%ld > MS[%ld]=%ld\n",
					(long)i, (long)pml[i], (long)i, (long)ms[i]);
			ret = 1;
		}
		if (pml[i] < 0) {
			fprintf(stderr, "  FAIL: PML[%ld]=%ld < 0\n", (long)i, (long)pml[i]);
			ret = 1;
		}
	}

	if (ret == 0) {
		fprintf(stderr, "  PML values:");
		for (i = 0; i < plen; ++i) fprintf(stderr, " %ld", (long)pml[i]);
		fprintf(stderr, "\n  PASS\n");
	}

	free(ms);
	free(pml);
	rb3_lcp_destroy(lcp);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Test: pattern matches entirely */
static int test_full_match(void)
{
	/* Text: ACGT, Pattern: ACG (should match fully with length 3 at pos 0) */
	int64_t text[] = {1, 2, 3, 4}; /* ACGT */
	uint8_t pat[] = {1, 2, 3}; /* ACG */
	return test_ms("full_match", text, 4, pat, 3);
}

/* Test: pattern not in text */
static int test_no_match(void)
{
	/* Text: AAAA, Pattern: CG (C and G don't appear in text) */
	int64_t text[] = {1, 1, 1, 1}; /* AAAA */
	uint8_t pat[] = {2, 3}; /* CG */
	return test_ms("no_match", text, 4, pat, 2);
}

/* Test: single char pattern */
static int test_single_char(void)
{
	int64_t text[] = {1, 2, 3, 4}; /* ACGT */
	uint8_t pat[] = {1}; /* A */
	return test_ms("single_char", text, 4, pat, 1);
}

/* Test: pattern with partial matches */
static int test_partial_match(void)
{
	/* Text: ACGTACGT, Pattern: ACGA (ACG matches but A after G doesn't continue as ACGT) */
	int64_t text[] = {1, 2, 3, 4, 1, 2, 3, 4}; /* ACGTACGT */
	uint8_t pat[] = {1, 2, 3, 1}; /* ACGA */
	return test_ms("partial_match", text, 8, pat, 4);
}

/* Test: repetitive text and pattern */
static int test_repetitive(void)
{
	/* Text: ACACAC, Pattern: ACAC */
	int64_t text[] = {1, 2, 1, 2, 1, 2}; /* ACACAC */
	uint8_t pat[] = {1, 2, 1, 2}; /* ACAC */
	return test_ms("repetitive", text, 6, pat, 4);
}

/* Test: longer pattern than text */
static int test_long_pattern(void)
{
	int64_t text[] = {1, 2}; /* AC */
	uint8_t pat[] = {1, 2, 1, 2, 1}; /* ACACA */
	return test_ms("long_pattern", text, 2, pat, 5);
}

/* Test: varied text with complex pattern */
static int test_complex(void)
{
	/* Text: ACGTTTAAACCCCGGGG */
	int64_t text[] = {1, 2, 3, 4, 4, 4, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
	/* Pattern: TTTAAGCC */
	uint8_t pat[] = {4, 4, 4, 1, 1, 3, 2, 2};
	return test_ms("complex", text, 17, pat, 8);
}

/* Test: PML properties */
static int test_pml_properties(void)
{
	int64_t text[] = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4}; /* ACGTACGTACGT */
	uint8_t pat[] = {1, 2, 3, 1, 4, 2, 3, 4}; /* ACGATCGT */
	int ret = 0;
	ret |= test_pml("pml_repetitive", text, 12, pat, 8);

	{
		int64_t text2[] = {1, 2, 1, 2, 1, 2}; /* ACACAC */
		uint8_t pat2[] = {1, 2, 1, 3}; /* ACAG */
		ret |= test_pml("pml_partial", text2, 6, pat2, 4);
	}
	return ret;
}

/* Test: empty and edge cases */
static int test_edge_cases(void)
{
	rb3_fmi_t fmi = {0};
	rb3_lcp_t *lcp;
	int64_t text[] = {1, 2, 3};
	int64_t ms[1];
	int ret = 0;

	fprintf(stderr, "Test MS: edge_cases\n");

	if (build_fmi_from_text(&fmi, text, 3) < 0) {
		fprintf(stderr, "  FAIL: could not build FMI\n");
		return 1;
	}
	lcp = rb3_lcp_build(&fmi);
	rb3_lcp_build_thresholds(lcp);

	/* Empty pattern should not crash */
	rb3_ms_compute(&fmi, lcp, (const uint8_t *)"", 0, ms);
	rb3_pml_compute(&fmi, lcp, (const uint8_t *)"", 0, ms);

	fprintf(stderr, "  PASS\n");

	rb3_lcp_destroy(lcp);
	rb3_fmi_free(&fmi);
	return ret;
}

int main(void)
{
	int ret = 0;

	rb3_verbose = 3;
	fprintf(stderr, "=== Matching Statistics Tests ===\n\n");

	ret |= test_edge_cases();
	ret |= test_single_char();
	ret |= test_full_match();
	ret |= test_no_match();
	ret |= test_partial_match();
	ret |= test_repetitive();
	ret |= test_long_pattern();
	ret |= test_complex();
	ret |= test_pml_properties();

	fprintf(stderr, "\n");
	if (ret == 0)
		fprintf(stderr, "=== ALL TESTS PASSED ===\n");
	else
		fprintf(stderr, "=== SOME TESTS FAILED ===\n");

	return ret;
}
