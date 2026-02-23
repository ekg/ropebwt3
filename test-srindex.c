#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "srindex.h"

extern int rb3_verbose;
extern int rb3_dbg_flag;

/* Build a plain BWT from a text using suffix array construction. */
static void build_sa_naive(const uint8_t *text, int64_t n, int64_t *sa)
{
	int64_t i, j;
	for (i = 0; i < n; ++i)
		sa[i] = i;
	for (i = 0; i < n; ++i) {
		for (j = i + 1; j < n; ++j) {
			int64_t a = sa[i], b = sa[j], k;
			for (k = 0; k < n; ++k) {
				if (text[(a + k) % n] < text[(b + k) % n]) break;
				if (text[(a + k) % n] > text[(b + k) % n]) {
					sa[i] = b;
					sa[j] = a;
					break;
				}
			}
		}
	}
}

static void build_bwt_from_sa(const uint8_t *text, int64_t n, const int64_t *sa, uint8_t *bwt)
{
	int64_t i;
	for (i = 0; i < n; ++i)
		bwt[i] = text[(sa[i] + n - 1) % n];
}

static int64_t count_runs(const uint8_t *bwt, int64_t n)
{
	int64_t i, r = 1;
	for (i = 1; i < n; ++i)
		if (bwt[i] != bwt[i-1]) ++r;
	return r;
}

/*
 * Test the SR-index with a given string and subsampling parameter s.
 * Verifies:
 * 1. phi(SA[k]) = SA[k-1] for all k in [1, n)
 * 2. Toehold correctness at run boundaries
 * 3. locate (original API with provided toehold)
 * 4. locate_one (LF-walk to subsampled sample)
 * 5. locate_all (automatic toehold resolution + phi)
 * 6. Space usage
 */
static int test_string(const char *name, const uint8_t *text, int64_t n, int32_t s)
{
	int64_t *sa, i;
	uint8_t *bwt;
	rld_t *e;
	rb3_fmi_t fmi;
	rb3_srindex_t *sr;
	int errors = 0;
	int64_t r;

	printf("=== %s (n=%lld, s=%d) ===\n", name, (long long)n, s);

	sa = (int64_t*)malloc(n * sizeof(int64_t));
	bwt = (uint8_t*)malloc(n);

	build_sa_naive(text, n, sa);
	build_bwt_from_sa(text, n, sa, bwt);
	r = count_runs(bwt, n);

	printf("n=%lld, r=%lld\n", (long long)n, (long long)r);

	e = rb3_enc_plain2rld(n, bwt, 3);
	if (e == 0) {
		fprintf(stderr, "Failed to build RLD\n");
		free(sa); free(bwt);
		return 1;
	}
	memset(&fmi, 0, sizeof(fmi));
	rb3_fmi_init(&fmi, e, 0);

	/* Build SR-index with parameter s */
	sr = rb3_srindex_build(&fmi, s, 1);
	if (sr == 0) {
		fprintf(stderr, "Failed to build SR-index\n");
		rb3_fmi_free(&fmi);
		free(sa); free(bwt);
		return 1;
	}

	printf("SR-index: n_runs=%lld, n_sub=%lld, s=%d\n",
	       (long long)sr->n_runs, (long long)sr->n_sub, sr->s);

	/* 1. Verify phi function */
	{
		int phi_errors = 0;
		for (i = 1; i < n; ++i) {
			int64_t phi_val = rb3_srindex_phi(sr, sa[i]);
			if (phi_val != sa[i-1]) {
				if (phi_errors < 5)
					fprintf(stderr, "  phi(SA[%lld]=%lld) = %lld, expected %lld\n",
					        (long long)i, (long long)sa[i], (long long)phi_val, (long long)sa[i-1]);
				phi_errors++;
			}
		}
		if (phi_errors) {
			fprintf(stderr, "FAILED: %d phi mismatches\n", phi_errors);
			errors++;
		} else {
			printf("Phi function: OK\n");
		}
	}

	/* 2. Verify toehold at run boundaries */
	{
		int th_errors = 0;
		for (i = 0; i < sr->n_samples; ++i) {
			int64_t pos = sr->run_pos[i];
			int64_t th = rb3_srindex_toehold(sr, pos);
			if (th != sa[pos]) {
				if (th_errors < 5)
					fprintf(stderr, "  toehold(%lld) = %lld, expected SA[%lld]=%lld\n",
					        (long long)pos, (long long)th, (long long)pos, (long long)sa[pos]);
				th_errors++;
			}
		}
		if (th_errors) {
			fprintf(stderr, "FAILED: %d toehold mismatches\n", th_errors);
			errors++;
		} else {
			printf("Toehold: OK\n");
		}
	}

	/* 3. Verify locate (original API) with known toehold */
	{
		int64_t *out = (int64_t*)malloc(n * sizeof(int64_t));
		int64_t toehold = sa[n - 1];
		int64_t cnt = rb3_srindex_locate(sr, 0, n, toehold, out);
		if (cnt == n) {
			int loc_errors = 0;
			for (i = 0; i < n; ++i)
				if (out[i] != sa[i]) loc_errors++;
			if (loc_errors) {
				fprintf(stderr, "FAILED: locate had %d mismatches\n", loc_errors);
				errors++;
			} else {
				printf("Locate (toehold API): OK\n");
			}
		} else {
			fprintf(stderr, "FAILED: locate returned %lld, expected %lld\n",
			        (long long)cnt, (long long)n);
			errors++;
		}
		free(out);
	}

	/* 4. Verify locate_one: resolve SA[k] for all run boundary positions */
	{
		int lo_errors = 0, lo_tested = 0;
		for (i = 0; i < sr->n_samples; ++i) {
			int64_t pos = sr->run_pos[i];
			int64_t result = rb3_srindex_locate_one(sr, &fmi, pos);
			if (result != sa[pos]) {
				if (lo_errors < 5)
					fprintf(stderr, "  locate_one(%lld) = %lld, expected SA[%lld]=%lld\n",
					        (long long)pos, (long long)result, (long long)pos, (long long)sa[pos]);
				lo_errors++;
			}
			lo_tested++;
		}
		if (lo_errors) {
			fprintf(stderr, "FAILED: locate_one had %d mismatches out of %d\n",
			        lo_errors, lo_tested);
			errors++;
		} else {
			printf("Locate_one (run boundaries): OK (%d tested)\n", lo_tested);
		}
	}

	/* 5. Verify locate_all: full interval [0, n) */
	{
		int64_t *out = (int64_t*)malloc(n * sizeof(int64_t));
		int64_t cnt = rb3_srindex_locate_all(sr, &fmi, 0, n, out, n);
		if (cnt == n) {
			int la_errors = 0;
			for (i = 0; i < n; ++i)
				if (out[i] != sa[i]) la_errors++;
			if (la_errors) {
				fprintf(stderr, "FAILED: locate_all had %d mismatches\n", la_errors);
				errors++;
			} else {
				printf("Locate_all [0,%lld): OK\n", (long long)n);
			}
		} else {
			fprintf(stderr, "FAILED: locate_all returned %lld, expected %lld\n",
			        (long long)cnt, (long long)n);
			errors++;
		}
		free(out);
	}

	/* 5b. Verify locate_all with sub-intervals */
	if (n >= 6) {
		int64_t lo2 = 2, hi2 = n < 8 ? n - 1 : 8;
		int64_t *out = (int64_t*)malloc(n * sizeof(int64_t));
		int64_t cnt = rb3_srindex_locate_all(sr, &fmi, lo2, hi2, out, n);
		if (cnt == hi2 - lo2) {
			int la_errors = 0;
			for (i = 0; i < cnt; ++i)
				if (out[i] != sa[lo2 + i]) la_errors++;
			if (la_errors) {
				fprintf(stderr, "FAILED: locate_all sub-interval had %d mismatches\n", la_errors);
				errors++;
			} else {
				printf("Locate_all [%lld,%lld): OK\n", (long long)lo2, (long long)hi2);
			}
		} else {
			fprintf(stderr, "FAILED: locate_all sub-interval returned %lld, expected %lld\n",
			        (long long)cnt, (long long)(hi2 - lo2));
			errors++;
		}
		free(out);
	}

	/* 5c. Verify locate_one from ALL BWT positions (not just run boundaries) */
	{
		int lo_errors = 0;
		for (i = 0; i < n; ++i) {
			int64_t result = rb3_srindex_locate_one(sr, &fmi, i);
			if (result != sa[i]) {
				if (lo_errors < 5)
					fprintf(stderr, "  locate_one(%lld) = %lld, expected SA[%lld]=%lld\n",
					        (long long)i, (long long)result, (long long)i, (long long)sa[i]);
				lo_errors++;
			}
		}
		if (lo_errors) {
			fprintf(stderr, "FAILED: locate_one (all positions) had %d mismatches out of %lld\n",
			        lo_errors, (long long)n);
			errors++;
		} else {
			printf("Locate_one (all %lld positions): OK\n", (long long)n);
		}
	}

	/* 5d. Verify locate_all with max_pos limiting.
	 * When max_pos < hi-lo, locate_all returns the LAST max_pos positions
	 * (from toehold backwards), i.e. SA[hi-max_pos], ..., SA[hi-1]. */
	if (n >= 4) {
		int64_t max_p = 3;
		int64_t *out = (int64_t*)malloc(max_p * sizeof(int64_t));
		int64_t cnt = rb3_srindex_locate_all(sr, &fmi, 0, n, out, max_p);
		if (cnt == max_p) {
			int la_errors = 0;
			int64_t offset = n - max_p; /* SA[n-max_p], ..., SA[n-1] */
			for (i = 0; i < cnt; ++i)
				if (out[i] != sa[offset + i]) la_errors++;
			if (la_errors) {
				fprintf(stderr, "FAILED: locate_all max_pos had %d mismatches\n", la_errors);
				errors++;
			} else {
				printf("Locate_all max_pos=%lld: OK\n", (long long)max_p);
			}
		} else {
			fprintf(stderr, "FAILED: locate_all max_pos returned %lld, expected %lld\n",
			        (long long)cnt, (long long)max_p);
			errors++;
		}
		free(out);
	}

	/* 6. Verify space: subsampled sample count */
	{
		int64_t expected_sub;
		if (s == 1) {
			expected_sub = sr->n_samples; /* all run boundary samples */
		} else {
			expected_sub = n / s + (n % s ? 1 : 0); /* ~n/s positions with SA%s==0 */
		}
		printf("Space: n_sub=%lld, expected~=%lld, n_runs=%lld, n/s=%lld\n",
		       (long long)sr->n_sub, (long long)expected_sub,
		       (long long)sr->n_runs, (long long)(n / s));

		if (s > 1 && sr->n_sub != expected_sub) {
			fprintf(stderr, "FAILED: n_sub=%lld, expected=%lld\n",
			        (long long)sr->n_sub, (long long)expected_sub);
			errors++;
		}
		if (s > 1) {
			/* Verify all stored samples have SA % s == 0 */
			int bad = 0;
			for (i = 0; i < sr->n_sub; ++i)
				if (sr->sub_sa[i] % s != 0) bad++;
			if (bad) {
				fprintf(stderr, "FAILED: %d subsampled entries with SA %% s != 0\n", bad);
				errors++;
			}
		}
	}

	/* 7. Verify serialization round-trip (dump + restore) */
	{
		const char *tmpfn = "/tmp/test_srindex_roundtrip.sri";
		rb3_srindex_t *sr2;
		int rt_errors = 0;

		if (rb3_srindex_dump(sr, tmpfn) != 0) {
			fprintf(stderr, "FAILED: dump returned error\n");
			errors++;
		} else {
			sr2 = rb3_srindex_restore(tmpfn);
			if (sr2 == 0) {
				fprintf(stderr, "FAILED: restore returned NULL\n");
				errors++;
			} else {
				/* Compare all fields */
				if (sr2->n != sr->n) { fprintf(stderr, "  roundtrip: n mismatch\n"); rt_errors++; }
				if (sr2->n_runs != sr->n_runs) { fprintf(stderr, "  roundtrip: n_runs mismatch\n"); rt_errors++; }
				if (sr2->n_samples != sr->n_samples) { fprintf(stderr, "  roundtrip: n_samples mismatch\n"); rt_errors++; }
				if (sr2->n_sub != sr->n_sub) { fprintf(stderr, "  roundtrip: n_sub mismatch\n"); rt_errors++; }
				if (sr2->s != sr->s) { fprintf(stderr, "  roundtrip: s mismatch\n"); rt_errors++; }
				if (sr2->m != sr->m) { fprintf(stderr, "  roundtrip: m mismatch\n"); rt_errors++; }

				for (i = 0; i < sr->n_runs && rt_errors < 5; ++i) {
					if (sr2->phi_sa[i] != sr->phi_sa[i]) {
						fprintf(stderr, "  roundtrip: phi_sa[%lld] %lld != %lld\n",
						        (long long)i, (long long)sr2->phi_sa[i], (long long)sr->phi_sa[i]);
						rt_errors++;
					}
					if (sr2->phi_da[i] != sr->phi_da[i]) {
						fprintf(stderr, "  roundtrip: phi_da[%lld] %lld != %lld\n",
						        (long long)i, (long long)sr2->phi_da[i], (long long)sr->phi_da[i]);
						rt_errors++;
					}
				}
				for (i = 0; i < sr->n_samples && rt_errors < 5; ++i) {
					if (sr2->run_pos[i] != sr->run_pos[i]) {
						fprintf(stderr, "  roundtrip: run_pos[%lld] %lld != %lld\n",
						        (long long)i, (long long)sr2->run_pos[i], (long long)sr->run_pos[i]);
						rt_errors++;
					}
					if (sr2->run_sa[i] != sr->run_sa[i]) {
						fprintf(stderr, "  roundtrip: run_sa[%lld] %lld != %lld\n",
						        (long long)i, (long long)sr2->run_sa[i], (long long)sr->run_sa[i]);
						rt_errors++;
					}
				}
				if (!sr->sub_is_alias) {
					for (i = 0; i < sr->n_sub && rt_errors < 5; ++i) {
						if (sr2->sub_pos[i] != sr->sub_pos[i]) {
							fprintf(stderr, "  roundtrip: sub_pos[%lld] %lld != %lld\n",
							        (long long)i, (long long)sr2->sub_pos[i], (long long)sr->sub_pos[i]);
							rt_errors++;
						}
						if (sr2->sub_sa[i] != sr->sub_sa[i]) {
							fprintf(stderr, "  roundtrip: sub_sa[%lld] %lld != %lld\n",
							        (long long)i, (long long)sr2->sub_sa[i], (long long)sr->sub_sa[i]);
							rt_errors++;
						}
					}
				}
				/* Verify phi function on restored index */
				for (i = 1; i < n && rt_errors < 5; ++i) {
					int64_t phi_val = rb3_srindex_phi(sr2, sa[i]);
					if (phi_val != sa[i-1]) {
						fprintf(stderr, "  roundtrip phi(SA[%lld]=%lld) = %lld, expected %lld\n",
						        (long long)i, (long long)sa[i], (long long)phi_val, (long long)sa[i-1]);
						rt_errors++;
					}
				}
				/* Verify locate_all on restored index */
				{
					int64_t *out2 = (int64_t*)malloc(n * sizeof(int64_t));
					int64_t cnt2 = rb3_srindex_locate_all(sr2, &fmi, 0, n, out2, n);
					if (cnt2 != n) {
						fprintf(stderr, "  roundtrip locate_all returned %lld, expected %lld\n",
						        (long long)cnt2, (long long)n);
						rt_errors++;
					} else {
						for (i = 0; i < n && rt_errors < 5; ++i)
							if (out2[i] != sa[i]) rt_errors++;
					}
					free(out2);
				}

				if (rt_errors) {
					fprintf(stderr, "FAILED: serialization roundtrip had %d mismatches\n", rt_errors);
					errors++;
				} else {
					printf("Serialization roundtrip: OK\n");
				}
				rb3_srindex_destroy(sr2);
			}
		}
		remove(tmpfn);
	}

	rb3_srindex_destroy(sr);
	rb3_fmi_free(&fmi);
	free(sa); free(bwt);

	if (errors == 0) {
		printf("PASSED\n\n");
		return 0;
	} else {
		printf("FAILED (%d errors)\n\n", errors);
		return 1;
	}
}

/* Test string definitions */
static int test_small(int32_t s)
{
	uint8_t text[] = {1, 1, 2, 3, 0}; /* AACG$ */
	return test_string("Small AACG$", text, 5, s);
}

static int test_repetitive(int32_t s)
{
	uint8_t text[] = {1,2,1,2,1,2,1,2,1,2,0}; /* (AC)^5$ */
	return test_string("Repetitive (AC)^5$", text, 11, s);
}

static int test_dna(int32_t s)
{
	uint8_t text[] = {1,2,3,4,1,2,3,4,0}; /* ACGTACGT$ */
	return test_string("DNA ACGTACGT$", text, 9, s);
}

static int test_longer(int32_t s)
{
	uint8_t text[] = {1,2,3,4,1,2,3,4,1,2,3,4,5,5,1,2,3,4,1,2,3,4,0};
	return test_string("Longer DNA with N's", text, 23, s);
}

/* Generate a longer repetitive string for testing large s values */
static int test_very_repetitive(int32_t s)
{
	/* "ACGTACGT" repeated 16 times + "$" = 129 chars */
	int64_t n = 129;
	uint8_t *text = (uint8_t*)malloc(n);
	int ret, i;
	uint8_t unit[] = {1,2,3,4,1,2,3,4};
	for (i = 0; i < 128; ++i)
		text[i] = unit[i % 8];
	text[128] = 0; /* sentinel */
	ret = test_string("Very repetitive (ACGTACGT)^16$", text, n, s);
	free(text);
	return ret;
}

/* Test with a string having many occurrences of a pattern */
static int test_many_occurrences(int32_t s)
{
	/* "AAAA...A$" with 100 A's + sentinel = 101 chars.
	 * Pattern "A" has 100 occurrences.
	 * Pattern "AA" has 99 occurrences, etc. */
	int64_t n = 101;
	uint8_t *text = (uint8_t*)malloc(n);
	int ret, i;
	for (i = 0; i < 100; ++i) text[i] = 1; /* A */
	text[100] = 0; /* sentinel */
	ret = test_string("Many A's (100)", text, n, s);
	free(text);
	return ret;
}

static int cmp_int64_t(const void *a, const void *b)
{
	int64_t x = *(const int64_t*)a;
	int64_t y = *(const int64_t*)b;
	return (x > y) - (x < y);
}

/*
 * Test pattern search + locate with 1, 10, 100, 1000 occurrences.
 * Uses text "A"x1000 + "$" (n=1001) whose BWT is [A x 1000, $] with 2 runs.
 * Pattern "A"xm has exactly (1001-m) occurrences at positions {0, 1, ..., 1000-m}.
 */
static int test_large_pattern_locate(int32_t s)
{
	int64_t n = 1001, i;
	uint8_t *bwt;
	rld_t *e;
	rb3_fmi_t fmi;
	rb3_srindex_t *sr;
	int errors = 0;
	int pat_lens[] = {1000, 991, 901, 1};
	int pat_occs[] = {1, 10, 100, 1000};
	int t;

	printf("=== Pattern locate (n=%lld, s=%d) ===\n", (long long)n, s);

	/* Build BWT directly: 1000 A's then 1 $ */
	bwt = (uint8_t*)malloc(n);
	for (i = 0; i < 1000; ++i) bwt[i] = 1;
	bwt[1000] = 0;

	e = rb3_enc_plain2rld(n, bwt, 3);
	free(bwt);
	if (e == 0) {
		fprintf(stderr, "Failed to build RLD\n");
		return 1;
	}
	memset(&fmi, 0, sizeof(fmi));
	rb3_fmi_init(&fmi, e, 0);

	sr = rb3_srindex_build(&fmi, s, 1);
	if (sr == 0) {
		fprintf(stderr, "Failed to build SR-index\n");
		rb3_fmi_free(&fmi);
		return 1;
	}

	printf("SR-index: n_runs=%lld, n_sub=%lld, s=%d\n",
	       (long long)sr->n_runs, (long long)sr->n_sub, sr->s);

	for (t = 0; t < 4; ++t) {
		int m = pat_lens[t];
		int expected = pat_occs[t];
		int64_t lo, hi, occ, cnt, *positions;
		int j;

		/* Backward search for "A" x m */
		lo = fmi.acc[1];
		hi = fmi.acc[2];
		for (j = 1; j < m; ++j)
			rb3_fmi_extend1(&fmi, &lo, &hi, 1);
		occ = hi - lo;

		if (occ != expected) {
			fprintf(stderr, "Pattern A*%d: expected %d occ, got %lld\n",
			        m, expected, (long long)occ);
			errors++;
			continue;
		}

		positions = (int64_t*)malloc(occ * sizeof(int64_t));
		cnt = rb3_srindex_locate_all(sr, &fmi, lo, hi, positions, occ);
		if (cnt != occ) {
			fprintf(stderr, "Pattern A*%d: locate_all returned %lld, expected %lld\n",
			        m, (long long)cnt, (long long)occ);
			errors++;
			free(positions);
			continue;
		}

		/* Sort and verify: expected positions are {0, 1, ..., 1000-m} */
		qsort(positions, cnt, sizeof(int64_t), cmp_int64_t);
		{
			int mismatch = 0;
			for (i = 0; i < cnt; ++i) {
				if (positions[i] != i) mismatch++;
			}
			if (mismatch) {
				fprintf(stderr, "Pattern A*%d: %d position mismatches\n", m, mismatch);
				errors++;
			} else {
				printf("Pattern A*%d: %d occ, locate OK\n", m, (int)occ);
			}
		}
		free(positions);
	}

	/* Space check */
	if (s > 1) {
		int64_t expected_sub = (n - 1) / s + 1;
		if (sr->n_sub != expected_sub) {
			fprintf(stderr, "Space: n_sub=%lld != expected %lld\n",
			        (long long)sr->n_sub, (long long)expected_sub);
			errors++;
		}
	}

	rb3_srindex_destroy(sr);
	rb3_fmi_free(&fmi);

	if (errors == 0) {
		printf("PASSED\n\n");
		return 0;
	} else {
		printf("FAILED (%d errors)\n\n", errors);
		return 1;
	}
}

int main(void)
{
	int ret = 0;
	int32_t s_values[] = {1, 4, 16, 64};
	int n_s = 4;
	int si;

	rb3_verbose = 3;

	for (si = 0; si < n_s; ++si) {
		int32_t s = s_values[si];
		printf("========== Testing with s=%d ==========\n\n", s);
		ret |= test_small(s);
		ret |= test_repetitive(s);
		ret |= test_dna(s);
		ret |= test_longer(s);
		ret |= test_very_repetitive(s);
		ret |= test_many_occurrences(s);
		ret |= test_large_pattern_locate(s);
	}

	if (ret == 0)
		printf("ALL TESTS PASSED\n");
	else
		printf("SOME TESTS FAILED\n");
	return ret;
}
