#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "move.h"

/* Helper: compute rank-based LF-mapping for any BWT position */
static int64_t rank_lf(const rb3_fmi_t *fmi, int64_t pos)
{
	int64_t ok[RB3_ASIZE];
	int c = rb3_fmi_rank1a(fmi, pos, ok);
	return fmi->acc[c] + ok[c];
}

/* Helper: find the run containing BWT position pos (linear scan) */
static int64_t find_run(const rb3_move_t *m, int64_t pos)
{
	int64_t i;
	for (i = 0; i < m->n_runs; i++)
		if (pos >= m->rows[i].p && pos < m->rows[i].p + m->rows[i].len)
			return i;
	return -1;
}

static int test_basic(void)
{
	/*
	 * Test BWT: [2, 1, 1, 0, 2, 1, 4, 4, 1, 2]
	 * Counts: {$:1, A:4, C:3, G:0, T:2, N:0}
	 * acc[]: [0, 1, 5, 8, 8, 10, 10]
	 *
	 * Maximal runs (8 runs):
	 *   i=0: c=2, len=1, p=0,  pi=5,  xi=4
	 *   i=1: c=1, len=2, p=1,  pi=1,  xi=1
	 *   i=2: c=0, len=1, p=3,  pi=0,  xi=0
	 *   i=3: c=2, len=1, p=4,  pi=6,  xi=5
	 *   i=4: c=1, len=1, p=5,  pi=3,  xi=2
	 *   i=5: c=4, len=2, p=6,  pi=8,  xi=6
	 *   i=6: c=1, len=1, p=8,  pi=4,  xi=3
	 *   i=7: c=2, len=1, p=9,  pi=7,  xi=5
	 */
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	int64_t exp_acc[] = {0, 1, 5, 8, 8, 10, 10};
	int n_exp = 8;
	int8_t  exp_c[]   = {2, 1, 0, 2, 1, 4, 1, 2};
	int64_t exp_len[] = {1, 2, 1, 1, 1, 2, 1, 1};
	int64_t exp_p[]   = {0, 1, 3, 4, 5, 6, 8, 9};
	int64_t exp_pi[]  = {5, 1, 0, 6, 3, 8, 4, 7};
	int64_t exp_xi[]  = {4, 1, 0, 5, 2, 6, 3, 5};

	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int i, ret = 0;

	/* Build FMD from plain BWT */
	e = rb3_enc_plain2rld(len, bwt, 3);
	if (e == 0) { fprintf(stderr, "FAIL: rb3_enc_plain2rld returned NULL\n"); return 1; }
	rb3_fmi_init(&fmi, e, 0);

	/* Build move table */
	m = rb3_move_build(&fmi);
	if (m == 0) { fprintf(stderr, "FAIL: rb3_move_build returned NULL\n"); ret = 1; goto done; }

	/* Verify BWT length */
	if (m->bwt_len != len) {
		fprintf(stderr, "FAIL: bwt_len = %ld, expected %ld\n", (long)m->bwt_len, (long)len);
		ret = 1; goto done;
	}

	/* Verify acc[] */
	for (i = 0; i <= RB3_ASIZE; ++i) {
		if (m->acc[i] != exp_acc[i]) {
			fprintf(stderr, "FAIL: acc[%d] = %ld, expected %ld\n", i, (long)m->acc[i], (long)exp_acc[i]);
			ret = 1; goto done;
		}
	}

	/* Verify number of runs */
	if (m->n_runs != n_exp) {
		fprintf(stderr, "FAIL: n_runs = %ld, expected %d\n", (long)m->n_runs, n_exp);
		ret = 1; goto done;
	}

	/* Verify each row's fields */
	for (i = 0; i < n_exp; ++i) {
		rb3_move_row_t *r = &m->rows[i];
		if (r->c != exp_c[i] || r->len != exp_len[i] || r->p != exp_p[i] ||
		    r->pi != exp_pi[i] || r->xi != exp_xi[i]) {
			fprintf(stderr, "FAIL: row[%d] = (c=%d, len=%ld, p=%ld, pi=%ld, xi=%ld)\n"
				"      expected (c=%d, len=%ld, p=%ld, pi=%ld, xi=%ld)\n",
				i, r->c, (long)r->len, (long)r->p, (long)r->pi, (long)r->xi,
				exp_c[i], (long)exp_len[i], (long)exp_p[i], (long)exp_pi[i], (long)exp_xi[i]);
			ret = 1; goto done;
		}
	}

	/* Verify LF-mapping of run heads against rb3_fmi_rank1a */
	for (i = 0; i < n_exp; ++i) {
		int64_t ok[RB3_ASIZE];
		int c_at_p = rb3_fmi_rank1a(&fmi, m->rows[i].p, ok);
		int64_t pi_check = fmi.acc[c_at_p] + ok[c_at_p];
		if (c_at_p != m->rows[i].c) {
			fprintf(stderr, "FAIL: rank1a at row[%d].p=%ld returned c=%d, expected %d\n",
				i, (long)m->rows[i].p, c_at_p, m->rows[i].c);
			ret = 1; goto done;
		}
		if (pi_check != m->rows[i].pi) {
			fprintf(stderr, "FAIL: rank1a-based pi for row[%d] = %ld, expected %ld\n",
				i, (long)pi_check, (long)m->rows[i].pi);
			ret = 1; goto done;
		}
	}

	/* Verify destination indices: pi must fall within the run at rows[xi] */
	for (i = 0; i < n_exp; ++i) {
		int64_t xi = m->rows[i].xi;
		int64_t pi = m->rows[i].pi;
		if (pi < m->rows[xi].p || pi >= m->rows[xi].p + m->rows[xi].len) {
			fprintf(stderr, "FAIL: row[%d].pi=%ld not in run rows[%ld] = [%ld, %ld)\n",
				i, (long)pi, (long)xi, (long)m->rows[xi].p,
				(long)(m->rows[xi].p + m->rows[xi].len));
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_basic: PASS (%d runs verified)\n", n_exp);
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

static int test_single_char(void)
{
	/* BWT with a single character: [1, 1, 1, 1] — 1 run */
	uint8_t bwt[] = {1, 1, 1, 1};
	int64_t len = 4;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);

	if (m->n_runs != 1) {
		fprintf(stderr, "FAIL: single_char n_runs = %ld, expected 1\n", (long)m->n_runs);
		ret = 1; goto done;
	}
	if (m->rows[0].c != 1 || m->rows[0].len != 4 || m->rows[0].p != 0) {
		fprintf(stderr, "FAIL: single_char row[0] mismatch\n");
		ret = 1; goto done;
	}
	if (m->rows[0].xi != 0) {
		fprintf(stderr, "FAIL: single_char xi = %ld, expected 0\n", (long)m->rows[0].xi);
		ret = 1; goto done;
	}
	fprintf(stderr, "test_single_char: PASS\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

static int test_alternating(void)
{
	/* Alternating BWT: [1, 2, 1, 2, 1, 2] — 6 runs, each length 1 */
	uint8_t bwt[] = {1, 2, 1, 2, 1, 2};
	int64_t len = 6;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int i, ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);

	if (m->n_runs != 6) {
		fprintf(stderr, "FAIL: alternating n_runs = %ld, expected 6\n", (long)m->n_runs);
		ret = 1; goto done;
	}

	/* Verify every run has length 1 */
	for (i = 0; i < 6; ++i) {
		if (m->rows[i].len != 1) {
			fprintf(stderr, "FAIL: alternating row[%d].len = %ld, expected 1\n",
				i, (long)m->rows[i].len);
			ret = 1; goto done;
		}
	}

	/* Verify destination indices */
	for (i = 0; i < 6; ++i) {
		int64_t xi = m->rows[i].xi;
		int64_t pi = m->rows[i].pi;
		if (pi < m->rows[xi].p || pi >= m->rows[xi].p + m->rows[xi].len) {
			fprintf(stderr, "FAIL: alternating row[%d].pi=%ld not in run at xi=%ld\n",
				i, (long)pi, (long)xi);
			ret = 1; goto done;
		}
	}

	/* Verify LF-mapping against rank queries */
	for (i = 0; i < 6; ++i) {
		int64_t ok[RB3_ASIZE];
		int c = rb3_fmi_rank1a(&fmi, m->rows[i].p, ok);
		int64_t pi = fmi.acc[c] + ok[c];
		if (pi != m->rows[i].pi) {
			fprintf(stderr, "FAIL: alternating rank check failed for row[%d]\n", i);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_alternating: PASS\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

static int test_fmr_backend(void)
{
	/* Test with FMR (rope) backend to verify both code paths */
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	mrope_t *r;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int i, ret = 0;

	r = rb3_enc_plain2fmr(len, bwt, 0, 0, 1);
	if (r == 0) { fprintf(stderr, "FAIL: rb3_enc_plain2fmr returned NULL\n"); return 1; }
	rb3_fmi_init(&fmi, 0, r);

	m = rb3_move_build(&fmi);
	if (m == 0) { fprintf(stderr, "FAIL: rb3_move_build (FMR) returned NULL\n"); ret = 1; goto done; }

	if (m->n_runs != 8) {
		fprintf(stderr, "FAIL: FMR n_runs = %ld, expected 8\n", (long)m->n_runs);
		ret = 1; goto done;
	}

	/* Verify LF-mapping and destination indices against rank queries */
	for (i = 0; i < m->n_runs; ++i) {
		int64_t ok[RB3_ASIZE];
		int c = rb3_fmi_rank1a(&fmi, m->rows[i].p, ok);
		int64_t pi = fmi.acc[c] + ok[c];
		int64_t xi = m->rows[i].xi;

		if (c != m->rows[i].c) {
			fprintf(stderr, "FAIL: FMR row[%d] c=%d, rank says %d\n", i, m->rows[i].c, c);
			ret = 1; goto done;
		}
		if (pi != m->rows[i].pi) {
			fprintf(stderr, "FAIL: FMR row[%d] pi=%ld, rank says %ld\n",
				i, (long)m->rows[i].pi, (long)pi);
			ret = 1; goto done;
		}
		if (m->rows[i].pi < m->rows[xi].p ||
		    m->rows[i].pi >= m->rows[xi].p + m->rows[xi].len) {
			fprintf(stderr, "FAIL: FMR row[%d] xi check failed\n", i);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_fmr_backend: PASS\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test LF-mapping via move table against rank-based LF for every position.
 */
static int test_lf_all_positions(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int64_t pos;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);

	for (pos = 0; pos < len; pos++) {
		int64_t run_idx = find_run(m, pos);
		int64_t move_lf = rb3_move_lf(m, pos, &run_idx);
		int64_t rank_lf_val = rank_lf(&fmi, pos);
		if (move_lf != rank_lf_val) {
			fprintf(stderr, "FAIL: lf_all pos=%ld: move_lf=%ld, rank_lf=%ld\n",
				(long)pos, (long)move_lf, (long)rank_lf_val);
			ret = 1; goto done;
		}
		/* Verify run_idx is correct after LF */
		if (run_idx < 0 || run_idx >= m->n_runs ||
		    move_lf < m->rows[run_idx].p ||
		    move_lf >= m->rows[run_idx].p + m->rows[run_idx].len) {
			fprintf(stderr, "FAIL: lf_all pos=%ld: returned run_idx=%ld doesn't contain lf=%ld\n",
				(long)pos, (long)run_idx, (long)move_lf);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_lf_all_positions: PASS (%ld positions)\n", (long)len);
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test that run splitting doesn't change LF-mapping results.
 */
static int test_split_preserves_lf(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int64_t pos;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);

	/* Apply splitting with d=2 */
	rb3_move_split(m, 2);
	if (m->d != 2) {
		fprintf(stderr, "FAIL: split d=%d, expected 2\n", m->d);
		ret = 1; goto done;
	}

	/* Verify LF-mapping still matches rank-based for every position */
	for (pos = 0; pos < len; pos++) {
		int64_t run_idx = find_run(m, pos);
		int64_t move_lf;
		int64_t rank_lf_val;
		if (run_idx < 0) {
			fprintf(stderr, "FAIL: split_lf pos=%ld: no run found\n", (long)pos);
			ret = 1; goto done;
		}
		move_lf = rb3_move_lf(m, pos, &run_idx);
		rank_lf_val = rank_lf(&fmi, pos);
		if (move_lf != rank_lf_val) {
			fprintf(stderr, "FAIL: split_lf pos=%ld: move=%ld, rank=%ld (n_runs=%ld)\n",
				(long)pos, (long)move_lf, (long)rank_lf_val, (long)m->n_runs);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_split_preserves_lf: PASS (d=2, %ld runs after split)\n", (long)m->n_runs);
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test run splitting with deeper depth (d=3).
 */
static int test_split_d3(void)
{
	/* Use a BWT with longer runs to actually trigger splitting */
	uint8_t bwt[] = {1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 0, 4,4,4,4};
	int64_t len = 20;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int64_t pos, orig_runs;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	orig_runs = m->n_runs;

	rb3_move_split(m, 3);

	/* Splitting should increase or keep n_runs */
	if (m->n_runs < orig_runs) {
		fprintf(stderr, "FAIL: split_d3 n_runs=%ld < orig=%ld\n",
			(long)m->n_runs, (long)orig_runs);
		ret = 1; goto done;
	}

	/* Verify LF for all positions */
	for (pos = 0; pos < len; pos++) {
		int64_t run_idx = find_run(m, pos);
		int64_t move_lf, rank_lf_val;
		if (run_idx < 0) {
			fprintf(stderr, "FAIL: split_d3 pos=%ld no run\n", (long)pos);
			ret = 1; goto done;
		}
		move_lf = rb3_move_lf(m, pos, &run_idx);
		rank_lf_val = rank_lf(&fmi, pos);
		if (move_lf != rank_lf_val) {
			fprintf(stderr, "FAIL: split_d3 pos=%ld: move=%ld, rank=%ld\n",
				(long)pos, (long)move_lf, (long)rank_lf_val);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_split_d3: PASS (orig=%ld, after=%ld runs)\n",
		(long)orig_runs, (long)m->n_runs);
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test reposition: verify that for each row and each character c,
 * the reposition distance points to a valid run of character c.
 */
static int test_reposition(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int64_t i;
	int c, ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	for (i = 0; i < m->n_runs; i++) {
		for (c = 0; c < RB3_ASIZE; c++) {
			/* Skip characters that don't appear */
			if (m->acc[c] == m->acc[c + 1]) continue;
			int64_t target = rb3_move_reposition(m, i, (int8_t)c);
			if (target < 0 || target >= m->n_runs) {
				fprintf(stderr, "FAIL: reposition row=%ld c=%d -> %ld (out of range)\n",
					(long)i, c, (long)target);
				ret = 1; goto done;
			}
			if (m->rows[target].c != c) {
				fprintf(stderr, "FAIL: reposition row=%ld c=%d -> row %ld has c=%d\n",
					(long)i, c, (long)target, m->rows[target].c);
				ret = 1; goto done;
			}
			/* Verify it's the nearest: no closer run of character c exists */
			{
				int64_t dist = target - i;
				int64_t j;
				if (dist > 0) {
					for (j = i + 1; j < target; j++)
						if (m->rows[j].c == c) {
							fprintf(stderr, "FAIL: reposition row=%ld c=%d -> %ld but row %ld also has c=%d\n",
								(long)i, c, (long)target, (long)j, c);
							ret = 1; goto done;
						}
				} else if (dist < 0) {
					for (j = i - 1; j > target; j--)
						if (m->rows[j].c == c) {
							fprintf(stderr, "FAIL: reposition row=%ld c=%d -> %ld but row %ld also has c=%d (bwd)\n",
								(long)i, c, (long)target, (long)j, c);
							ret = 1; goto done;
						}
				}
			}
		}
	}

	fprintf(stderr, "test_reposition: PASS\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test move_step: verify that a full backward search chain using move_step
 * produces the same results as rank-based backward extension.
 */
static int test_move_step(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int64_t pos;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	/* For each starting position and each character c that exists,
	 * verify move_step produces valid LF-mapping */
	for (pos = 0; pos < len; pos++) {
		int c;
		for (c = 0; c < RB3_ASIZE; c++) {
			int64_t run_idx, new_pos;
			if (m->acc[c] == m->acc[c + 1]) continue; /* skip absent chars */
			run_idx = find_run(m, pos);
			new_pos = rb3_move_step(m, pos, &run_idx, (int8_t)c);
			/* new_pos should be a valid position */
			if (new_pos < 0 || new_pos >= len) {
				fprintf(stderr, "FAIL: move_step pos=%ld c=%d -> %ld (out of range)\n",
					(long)pos, c, (long)new_pos);
				ret = 1; goto done;
			}
			/* run_idx should contain new_pos */
			if (run_idx < 0 || run_idx >= m->n_runs ||
			    new_pos < m->rows[run_idx].p ||
			    new_pos >= m->rows[run_idx].p + m->rows[run_idx].len) {
				fprintf(stderr, "FAIL: move_step pos=%ld c=%d -> run_idx=%ld doesn't contain %ld\n",
					(long)pos, c, (long)run_idx, (long)new_pos);
				ret = 1; goto done;
			}
			/* When c matches the current BWT character, move_step should equal rank_lf */
			if (bwt[pos] == c) {
				int64_t expected = rank_lf(&fmi, pos);
				if (new_pos != expected) {
					fprintf(stderr, "FAIL: move_step pos=%ld c=%d (match): got %ld, expected %ld\n",
						(long)pos, c, (long)new_pos, (long)expected);
					ret = 1; goto done;
				}
			}
		}
	}

	fprintf(stderr, "test_move_step: PASS\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test move_step with splitting enabled.
 */
static int test_move_step_split(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int64_t pos;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_split(m, 2);
	rb3_move_precompute_dist(m);

	for (pos = 0; pos < len; pos++) {
		int c;
		for (c = 0; c < RB3_ASIZE; c++) {
			int64_t run_idx, new_pos;
			if (m->acc[c] == m->acc[c + 1]) continue;
			run_idx = find_run(m, pos);
			new_pos = rb3_move_step(m, pos, &run_idx, (int8_t)c);
			if (new_pos < 0 || new_pos >= len) {
				fprintf(stderr, "FAIL: step_split pos=%ld c=%d -> %ld\n",
					(long)pos, c, (long)new_pos);
				ret = 1; goto done;
			}
			if (run_idx < 0 || run_idx >= m->n_runs ||
			    new_pos < m->rows[run_idx].p ||
			    new_pos >= m->rows[run_idx].p + m->rows[run_idx].len) {
				fprintf(stderr, "FAIL: step_split run_idx invalid\n");
				ret = 1; goto done;
			}
			if (bwt[pos] == c) {
				int64_t expected = rank_lf(&fmi, pos);
				if (new_pos != expected) {
					fprintf(stderr, "FAIL: step_split pos=%ld c=%d: got %ld, expected %ld\n",
						(long)pos, c, (long)new_pos, (long)expected);
					ret = 1; goto done;
				}
			}
		}
	}

	fprintf(stderr, "test_move_step_split: PASS\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test chained backward search: apply move_step repeatedly to simulate
 * backward search for a pattern, and verify against rank-based search.
 */
static int test_backward_search_chain(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	/* Do 20 chained LF-mapping steps from position 0, tracking via both methods */
	{
		int64_t pos = 0, run_idx;
		int step;
		run_idx = find_run(m, pos);
		for (step = 0; step < 20; step++) {
			int64_t rank_pos = rank_lf(&fmi, pos);
			int64_t move_pos = rb3_move_lf(m, pos, &run_idx);
			if (move_pos != rank_pos) {
				fprintf(stderr, "FAIL: chain step %d: move=%ld, rank=%ld\n",
					step, (long)move_pos, (long)rank_pos);
				ret = 1; goto done;
			}
			pos = move_pos;
		}
	}

	fprintf(stderr, "test_backward_search_chain: PASS (20 steps)\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test save/load round-trip: build, save, load, verify all rows match,
 * and verify LF-mapping works on the loaded (mmap'd) index.
 */
static int test_save_load(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m, *m2;
	int64_t pos;
	int i, ret = 0;
	const char *tmpfn = "/tmp/test-move.mvi";

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	if (rb3_move_save(m, tmpfn) != 0) {
		fprintf(stderr, "FAIL: save_load save returned error\n");
		ret = 1; goto done;
	}

	m2 = rb3_move_load(tmpfn);
	if (m2 == 0) {
		fprintf(stderr, "FAIL: save_load load returned NULL\n");
		ret = 1; goto done;
	}

	/* Compare metadata */
	if (m2->n_runs != m->n_runs || m2->bwt_len != m->bwt_len || m2->d != m->d) {
		fprintf(stderr, "FAIL: save_load metadata mismatch\n");
		ret = 1; goto done2;
	}
	for (i = 0; i <= RB3_ASIZE; i++) {
		if (m2->acc[i] != m->acc[i]) {
			fprintf(stderr, "FAIL: save_load acc[%d] mismatch: %ld vs %ld\n",
				i, (long)m2->acc[i], (long)m->acc[i]);
			ret = 1; goto done2;
		}
	}

	/* Compare every row */
	for (i = 0; i < m->n_runs; i++) {
		rb3_move_row_t *a = &m->rows[i], *b = &m2->rows[i];
		if (a->c != b->c || a->len != b->len || a->p != b->p ||
		    a->pi != b->pi || a->xi != b->xi) {
			fprintf(stderr, "FAIL: save_load row[%d] field mismatch\n", i);
			ret = 1; goto done2;
		}
		if (memcmp(a->dist, b->dist, sizeof(a->dist)) != 0) {
			fprintf(stderr, "FAIL: save_load row[%d] dist mismatch\n", i);
			ret = 1; goto done2;
		}
	}

	/* Verify LF-mapping works on the loaded (mmap'd) index */
	for (pos = 0; pos < len; pos++) {
		int64_t run_idx = find_run(m2, pos);
		int64_t move_lf = rb3_move_lf(m2, pos, &run_idx);
		int64_t rank_lf_val = rank_lf(&fmi, pos);
		if (move_lf != rank_lf_val) {
			fprintf(stderr, "FAIL: save_load LF mismatch at pos %ld: %ld vs %ld\n",
				(long)pos, (long)move_lf, (long)rank_lf_val);
			ret = 1; goto done2;
		}
	}

	fprintf(stderr, "test_save_load: PASS\n");
done2:
	rb3_move_destroy(m2);
done:
	unlink(tmpfn);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test save/load with run splitting: verify the split index round-trips correctly
 * and LF-mapping still works on the loaded version.
 */
static int test_save_load_split(void)
{
	uint8_t bwt[] = {1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 0, 4,4,4,4};
	int64_t len = 20;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m, *m2;
	int64_t pos;
	int i, ret = 0;
	const char *tmpfn = "/tmp/test-move-split.mvi";

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_split(m, 2);
	rb3_move_precompute_dist(m);

	if (rb3_move_save(m, tmpfn) != 0) {
		fprintf(stderr, "FAIL: save_load_split save returned error\n");
		ret = 1; goto done;
	}

	m2 = rb3_move_load(tmpfn);
	if (m2 == 0) {
		fprintf(stderr, "FAIL: save_load_split load returned NULL\n");
		ret = 1; goto done;
	}

	if (m2->n_runs != m->n_runs || m2->d != m->d) {
		fprintf(stderr, "FAIL: save_load_split metadata mismatch (n_runs=%ld vs %ld, d=%d vs %d)\n",
			(long)m2->n_runs, (long)m->n_runs, m2->d, m->d);
		ret = 1; goto done2;
	}

	/* Compare all rows byte-by-byte */
	for (i = 0; i < m->n_runs; i++) {
		if (memcmp(&m->rows[i], &m2->rows[i], sizeof(rb3_move_row_t)) != 0) {
			fprintf(stderr, "FAIL: save_load_split row[%d] mismatch\n", i);
			ret = 1; goto done2;
		}
	}

	/* Verify LF on loaded index */
	for (pos = 0; pos < len; pos++) {
		int64_t run_idx = find_run(m2, pos);
		int64_t move_lf, rank_lf_val;
		if (run_idx < 0) {
			fprintf(stderr, "FAIL: save_load_split pos=%ld: no run found\n", (long)pos);
			ret = 1; goto done2;
		}
		move_lf = rb3_move_lf(m2, pos, &run_idx);
		rank_lf_val = rank_lf(&fmi, pos);
		if (move_lf != rank_lf_val) {
			fprintf(stderr, "FAIL: save_load_split LF mismatch at pos %ld\n", (long)pos);
			ret = 1; goto done2;
		}
	}

	fprintf(stderr, "test_save_load_split: PASS\n");
done2:
	rb3_move_destroy(m2);
done:
	unlink(tmpfn);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/* Helper: rank-based backward search count (reference implementation) */
static int64_t rank_count(const rb3_fmi_t *f, int len, const uint8_t *pattern)
{
	int64_t lo, hi;
	int i;
	if (len <= 0) return f->acc[RB3_ASIZE];
	lo = f->acc[pattern[len - 1]];
	hi = f->acc[pattern[len - 1] + 1];
	for (i = len - 2; i >= 0 && lo < hi; --i) {
		int c = pattern[i];
		int64_t ok[RB3_ASIZE], ol[RB3_ASIZE];
		rb3_fmi_rank2a(f, lo, hi, ok, ol);
		lo = f->acc[c] + ok[c];
		hi = f->acc[c] + ol[c];
	}
	return hi > lo ? hi - lo : 0;
}

/*
 * Test move_count against rank-based count for the test BWT.
 * BWT: [2, 1, 1, 0, 2, 1, 4, 4, 1, 2]
 * This encodes the text: T = $AACCACAAT (with sentinel)
 * Actually let's just test all patterns up to length 4.
 */
static int test_count_basic(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	/* Test all single-character patterns */
	{
		int c;
		for (c = 0; c < RB3_ASIZE; c++) {
			uint8_t pat[1] = {(uint8_t)c};
			int64_t mc = rb3_move_count(m, 1, pat);
			int64_t rc = rank_count(&fmi, 1, pat);
			if (mc != rc) {
				fprintf(stderr, "FAIL: count_basic single c=%d: move=%ld, rank=%ld\n",
					c, (long)mc, (long)rc);
				ret = 1; goto done;
			}
		}
	}

	/* Test all 2-character patterns */
	{
		int c1, c2;
		for (c1 = 0; c1 < RB3_ASIZE; c1++) {
			for (c2 = 0; c2 < RB3_ASIZE; c2++) {
				uint8_t pat[2] = {(uint8_t)c1, (uint8_t)c2};
				int64_t mc = rb3_move_count(m, 2, pat);
				int64_t rc = rank_count(&fmi, 2, pat);
				if (mc != rc) {
					fprintf(stderr, "FAIL: count_basic 2-char [%d,%d]: move=%ld, rank=%ld\n",
						c1, c2, (long)mc, (long)rc);
					ret = 1; goto done;
				}
			}
		}
	}

	/* Test all 3-character patterns */
	{
		int c1, c2, c3;
		for (c1 = 0; c1 < RB3_ASIZE; c1++)
			for (c2 = 0; c2 < RB3_ASIZE; c2++)
				for (c3 = 0; c3 < RB3_ASIZE; c3++) {
					uint8_t pat[3] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3};
					int64_t mc = rb3_move_count(m, 3, pat);
					int64_t rc = rank_count(&fmi, 3, pat);
					if (mc != rc) {
						fprintf(stderr, "FAIL: count_basic 3-char [%d,%d,%d]: move=%ld, rank=%ld\n",
							c1, c2, c3, (long)mc, (long)rc);
						ret = 1; goto done;
					}
				}
	}

	/* Test empty pattern */
	{
		int64_t mc = rb3_move_count(m, 0, NULL);
		int64_t rc = rank_count(&fmi, 0, NULL);
		if (mc != rc) {
			fprintf(stderr, "FAIL: count_basic empty: move=%ld, rank=%ld\n",
				(long)mc, (long)rc);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_count_basic: PASS (1+36+216 patterns, plus empty)\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test move_count with a larger, more realistic BWT.
 */
static int test_count_larger(void)
{
	/* BWT of "ACAACG$" (constructed manually):
	 * Actually, let's use a longer example for thorough testing.
	 * BWT: [1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 0, 4,4,4,4]
	 */
	uint8_t bwt[] = {1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 0, 4,4,4,4};
	int64_t len = 20;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	/* Test all 2-character patterns */
	{
		int c1, c2;
		for (c1 = 0; c1 < RB3_ASIZE; c1++)
			for (c2 = 0; c2 < RB3_ASIZE; c2++) {
				uint8_t pat[2] = {(uint8_t)c1, (uint8_t)c2};
				int64_t mc = rb3_move_count(m, 2, pat);
				int64_t rc = rank_count(&fmi, 2, pat);
				if (mc != rc) {
					fprintf(stderr, "FAIL: count_larger 2-char [%d,%d]: move=%ld, rank=%ld\n",
						c1, c2, (long)mc, (long)rc);
					ret = 1; goto done;
				}
			}
	}

	/* Test all 4-character patterns (6^4 = 1296) */
	{
		int c1, c2, c3, c4;
		int n_tested = 0;
		for (c1 = 0; c1 < RB3_ASIZE; c1++)
			for (c2 = 0; c2 < RB3_ASIZE; c2++)
				for (c3 = 0; c3 < RB3_ASIZE; c3++)
					for (c4 = 0; c4 < RB3_ASIZE; c4++) {
						uint8_t pat[4] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3, (uint8_t)c4};
						int64_t mc = rb3_move_count(m, 4, pat);
						int64_t rc = rank_count(&fmi, 4, pat);
						if (mc != rc) {
							fprintf(stderr, "FAIL: count_larger 4-char [%d,%d,%d,%d]: move=%ld, rank=%ld\n",
								c1, c2, c3, c4, (long)mc, (long)rc);
							ret = 1; goto done;
						}
						n_tested++;
					}
		fprintf(stderr, "test_count_larger: PASS (%d patterns up to length 4)\n", 36 + n_tested);
	}

done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test move_count with splitting enabled.
 * After splitting, the cumulative rank table should still produce correct results.
 */
static int test_count_split(void)
{
	uint8_t bwt[] = {1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 0, 4,4,4,4};
	int64_t len = 20;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_split(m, 2);
	rb3_move_precompute_dist(m);

	/* Test all 3-character patterns (6^3 = 216) */
	{
		int c1, c2, c3;
		for (c1 = 0; c1 < RB3_ASIZE; c1++)
			for (c2 = 0; c2 < RB3_ASIZE; c2++)
				for (c3 = 0; c3 < RB3_ASIZE; c3++) {
					uint8_t pat[3] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3};
					int64_t mc = rb3_move_count(m, 3, pat);
					int64_t rc = rank_count(&fmi, 3, pat);
					if (mc != rc) {
						fprintf(stderr, "FAIL: count_split [%d,%d,%d]: move=%ld, rank=%ld\n",
							c1, c2, c3, (long)mc, (long)rc);
						ret = 1; goto done;
					}
				}
	}

	fprintf(stderr, "test_count_split: PASS (d=2, %ld runs, 216 patterns)\n", (long)m->n_runs);
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test rb3_fmi_t integration: verify the mv field is properly initialized.
 */
static int test_fmi_mv_field(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);

	/* mv should be NULL after init */
	if (fmi.mv != 0) {
		fprintf(stderr, "FAIL: fmi_mv_field: mv not NULL after init\n");
		ret = 1; goto done;
	}

	/* Build and attach move structure */
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);
	fmi.mv = m;

	/* Verify we can use it for counting */
	{
		uint8_t pat[] = {1}; /* single 'A' */
		int64_t mc = rb3_move_count(fmi.mv, 1, pat);
		int64_t rc = rank_count(&fmi, 1, pat);
		if (mc != rc) {
			fprintf(stderr, "FAIL: fmi_mv_field: count mismatch: move=%ld, rank=%ld\n",
				(long)mc, (long)rc);
			ret = 1; goto done;
		}
	}

	fprintf(stderr, "test_fmi_mv_field: PASS\n");
done:
	fmi.mv = 0; /* detach before free */
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test b-move: verify rb3_bmove_extend produces identical results to rb3_fmd_extend
 * for all possible intervals on a small BWT.
 */
static int test_bmove_extend(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_bmove_t *bm;
	int ret = 0;
	int64_t lo, hi;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	bm = rb3_bmove_init(m);
	if (bm == 0) { fprintf(stderr, "FAIL: bmove_extend: rb3_bmove_init returned NULL\n"); ret = 1; goto done; }

	/* Test all intervals [lo, hi) where 0 <= lo < hi <= len, both directions */
	for (lo = 0; lo < len; lo++) {
		for (hi = lo + 1; hi <= len; hi++) {
			int dir;
			for (dir = 0; dir <= 1; dir++) {
				rb3_sai_t ik, fmd_ok[6], bm_ok[6];
				int c;
				ik.x[0] = lo;
				ik.x[1] = hi > len/2 ? lo : hi; /* arbitrary but consistent x[1] */
				ik.size = hi - lo;
				ik.info = 0;
				rb3_fmd_extend(&fmi, &ik, fmd_ok, dir);
				rb3_bmove_extend(bm, &ik, bm_ok, dir);
				for (c = 0; c < RB3_ASIZE; c++) {
					if (bm_ok[c].x[0] != fmd_ok[c].x[0] || bm_ok[c].x[1] != fmd_ok[c].x[1] || bm_ok[c].size != fmd_ok[c].size) {
						fprintf(stderr, "FAIL: bmove_extend [%ld,%ld) dir=%d c=%d: bm=(%ld,%ld,%ld) fmd=(%ld,%ld,%ld)\n",
							(long)lo, (long)hi, dir, c,
							(long)bm_ok[c].x[0], (long)bm_ok[c].x[1], (long)bm_ok[c].size,
							(long)fmd_ok[c].x[0], (long)fmd_ok[c].x[1], (long)fmd_ok[c].size);
						ret = 1; goto done;
					}
				}
			}
		}
	}

	fprintf(stderr, "test_bmove_extend: PASS (all intervals, both directions)\n");
done:
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test b-move extend with run splitting (verify cumrank still correct after split).
 */
static int test_bmove_extend_split(void)
{
	uint8_t bwt[] = {1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 0, 4,4,4,4};
	int64_t len = 20;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_bmove_t *bm;
	int ret = 0;
	int64_t lo, hi;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_split(m, 2);
	bm = rb3_bmove_init(m);
	if (bm == 0) { fprintf(stderr, "FAIL: bmove_extend_split: init returned NULL\n"); ret = 1; goto done; }

	for (lo = 0; lo < len; lo++) {
		for (hi = lo + 1; hi <= len; hi++) {
			int dir;
			for (dir = 0; dir <= 1; dir++) {
				rb3_sai_t ik, fmd_ok[6], bm_ok[6];
				int c;
				ik.x[0] = lo;
				ik.x[1] = lo;
				ik.size = hi - lo;
				ik.info = 0;
				rb3_fmd_extend(&fmi, &ik, fmd_ok, dir);
				rb3_bmove_extend(bm, &ik, bm_ok, dir);
				for (c = 0; c < RB3_ASIZE; c++) {
					if (bm_ok[c].x[0] != fmd_ok[c].x[0] || bm_ok[c].x[1] != fmd_ok[c].x[1] || bm_ok[c].size != fmd_ok[c].size) {
						fprintf(stderr, "FAIL: bmove_extend_split [%ld,%ld) dir=%d c=%d\n",
							(long)lo, (long)hi, dir, c);
						ret = 1; goto done;
					}
				}
			}
		}
	}

	fprintf(stderr, "test_bmove_extend_split: PASS (d=2, %ld runs)\n", (long)m->n_runs);
done:
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test b-move SMEM on a symmetric BWT, comparing against FMD SMEM.
 *
 * Symmetric BWT for text ACAC$GTGT$ (sequences ACAC and GTGT = rc(ACAC)):
 *   BWT = [4, 2, 2, 0, 1, 1, 4, 0, 3, 3]
 *   acc = [0, 2, 4, 6, 8, 10, 10]
 */
static int test_bmove_smem(void)
{
	uint8_t bwt[] = {4, 2, 2, 0, 1, 1, 4, 0, 3, 3};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m = 0;
	rb3_bmove_t *bm = 0;
	rb3_sai_v fmd_mem = {0,0,0}, bm_mem = {0,0,0};
	int ret = 0;

	/* Test queries: various patterns encoded as nt6 */
	uint8_t q1[] = {1,2,1,2}; /* ACAC */
	uint8_t q2[] = {4,3,4,3}; /* TGTG - rc(CACA) */
	uint8_t q3[] = {1,2};     /* AC */
	uint8_t q4[] = {1};       /* A */
	uint8_t q5[] = {1,2,1,2, 3,4,3,4}; /* ACACGTGT - longer query */
	struct { uint8_t *q; int len; const char *name; } queries[] = {
		{q1, 4, "ACAC"}, {q2, 4, "TGTG"}, {q3, 2, "AC"}, {q4, 1, "A"}, {q5, 8, "ACACGTGT"},
	};
	int nq = 5, qi;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);

	/* Verify symmetry */
	if (!rb3_fmi_is_symmetric(&fmi)) {
		fprintf(stderr, "FAIL: bmove_smem: BWT not symmetric\n");
		ret = 1; goto done;
	}

	m = rb3_move_build(&fmi);
	bm = rb3_bmove_init(m);

	for (qi = 0; qi < nq; qi++) {
		int64_t fmd_n, bm_n;
		int64_t j;

		/* Test both original and TG algorithms */
		/* Original SMEM */
		fmd_mem.n = 0;
		bm_mem.n = 0;
		fmd_n = rb3_fmd_smem(0, &fmi, queries[qi].len, queries[qi].q, &fmd_mem, 1, 1);
		bm_n = rb3_bmove_smem(0, bm, queries[qi].len, queries[qi].q, &bm_mem, 1, 1);
		if (fmd_n != bm_n) {
			fprintf(stderr, "FAIL: bmove_smem ORI query=%s: fmd_n=%ld, bm_n=%ld\n",
				queries[qi].name, (long)fmd_n, (long)bm_n);
			ret = 1; goto done;
		}
		for (j = 0; j < fmd_n; j++) {
			if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
				fprintf(stderr, "FAIL: bmove_smem ORI query=%s mem[%ld]: fmd=(%ld,%ld,0x%lx) bm=(%ld,%ld,0x%lx)\n",
					queries[qi].name, (long)j,
					(long)fmd_mem.a[j].x[0], (long)fmd_mem.a[j].size, (unsigned long)fmd_mem.a[j].info,
					(long)bm_mem.a[j].x[0], (long)bm_mem.a[j].size, (unsigned long)bm_mem.a[j].info);
				ret = 1; goto done;
			}
		}

		/* TG SMEM */
		fmd_mem.n = 0;
		bm_mem.n = 0;
		fmd_n = rb3_fmd_smem_TG(0, &fmi, queries[qi].len, queries[qi].q, &fmd_mem, 1, 1);
		bm_n = rb3_bmove_smem_TG(0, bm, queries[qi].len, queries[qi].q, &bm_mem, 1, 1);
		if (fmd_n != bm_n) {
			fprintf(stderr, "FAIL: bmove_smem TG query=%s: fmd_n=%ld, bm_n=%ld\n",
				queries[qi].name, (long)fmd_n, (long)bm_n);
			ret = 1; goto done;
		}
		for (j = 0; j < fmd_n; j++) {
			if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
				fprintf(stderr, "FAIL: bmove_smem TG query=%s mem[%ld] mismatch\n",
					queries[qi].name, (long)j);
				ret = 1; goto done;
			}
		}
	}

	fprintf(stderr, "test_bmove_smem: PASS (%d queries, ORI+TG)\n", nq);
done:
	free(fmd_mem.a);
	free(bm_mem.a);
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test b-move SMEM with a larger symmetric BWT and various min_occ/min_len.
 *
 * Symmetric BWT for text ACA$TGT$ (sequences ACA and TGT = rc(ACA)):
 *   BWT = [4, 1, 2, 0, 1, 4, 3, 0]
 *   acc = [0, 2, 4, 5, 6, 8, 8]
 */
static int test_bmove_smem_params(void)
{
	uint8_t bwt[] = {4, 1, 2, 0, 1, 4, 3, 0};
	int64_t len = 8;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m = 0;
	rb3_bmove_t *bm = 0;
	rb3_sai_v fmd_mem = {0,0,0}, bm_mem = {0,0,0};
	int ret = 0;

	/* Query patterns */
	uint8_t q1[] = {1,2,1}; /* ACA */
	uint8_t q2[] = {4,3,4}; /* TGT */
	uint8_t q3[] = {1,2,1,4,3,4}; /* ACATGT */
	struct { uint8_t *q; int len; } queries[] = {
		{q1, 3}, {q2, 3}, {q3, 6},
	};
	int nq = 3;
	int64_t min_occs[] = {1, 2};
	int64_t min_lens[] = {1, 2, 3};
	int qi, oi, li;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	if (!rb3_fmi_is_symmetric(&fmi)) {
		fprintf(stderr, "FAIL: bmove_smem_params: BWT not symmetric\n");
		ret = 1; goto done;
	}
	m = rb3_move_build(&fmi);
	bm = rb3_bmove_init(m);

	for (qi = 0; qi < nq; qi++) {
		for (oi = 0; oi < 2; oi++) {
			for (li = 0; li < 3; li++) {
				int64_t fmd_n, bm_n;
				int64_t j, min_occ = min_occs[oi], min_len = min_lens[li];

				/* Test TG algorithm */
				fmd_mem.n = 0;
				bm_mem.n = 0;
				fmd_n = rb3_fmd_smem_TG(0, &fmi, queries[qi].len, queries[qi].q, &fmd_mem, min_occ, min_len);
				bm_n = rb3_bmove_smem_TG(0, bm, queries[qi].len, queries[qi].q, &bm_mem, min_occ, min_len);
				if (fmd_n != bm_n) {
					fprintf(stderr, "FAIL: bmove_smem_params TG qi=%d occ=%ld len=%ld: fmd=%ld bm=%ld\n",
						qi, (long)min_occ, (long)min_len, (long)fmd_n, (long)bm_n);
					ret = 1; goto done;
				}
				for (j = 0; j < fmd_n; j++) {
					if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
						fprintf(stderr, "FAIL: bmove_smem_params TG qi=%d occ=%ld len=%ld mem[%ld] mismatch\n",
							qi, (long)min_occ, (long)min_len, (long)j);
						ret = 1; goto done;
					}
				}

				/* Test ORI algorithm */
				fmd_mem.n = 0;
				bm_mem.n = 0;
				fmd_n = rb3_fmd_smem(0, &fmi, queries[qi].len, queries[qi].q, &fmd_mem, min_occ, min_len);
				bm_n = rb3_bmove_smem(0, bm, queries[qi].len, queries[qi].q, &bm_mem, min_occ, min_len);
				if (fmd_n != bm_n) {
					fprintf(stderr, "FAIL: bmove_smem_params ORI qi=%d occ=%ld len=%ld: fmd=%ld bm=%ld\n",
						qi, (long)min_occ, (long)min_len, (long)fmd_n, (long)bm_n);
					ret = 1; goto done;
				}
				for (j = 0; j < fmd_n; j++) {
					if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
						fprintf(stderr, "FAIL: bmove_smem_params ORI qi=%d occ=%ld len=%ld mem[%ld] mismatch\n",
							qi, (long)min_occ, (long)min_len, (long)j);
						ret = 1; goto done;
					}
				}
			}
		}
	}

	fprintf(stderr, "test_bmove_smem_params: PASS (%d queries × 2 min_occ × 3 min_len, ORI+TG)\n", nq);
done:
	free(fmd_mem.a);
	free(bm_mem.a);
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test b-move SMEM with exhaustive short queries on symmetric BWT.
 * Generate all possible 2-4 character queries (using nt6 chars 1-4)
 * and verify bmove matches fmd for both algorithms.
 */
static int test_bmove_smem_exhaustive(void)
{
	uint8_t bwt[] = {4, 2, 2, 0, 1, 1, 4, 0, 3, 3};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_bmove_t *bm;
	rb3_sai_v fmd_mem = {0,0,0}, bm_mem = {0,0,0};
	int ret = 0, n_tested = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	bm = rb3_bmove_init(m);

	/* All 2-char queries using chars 1-4 (A,C,G,T): 4^2 = 16 */
	{
		int c1, c2;
		for (c1 = 1; c1 <= 4; c1++) {
			for (c2 = 1; c2 <= 4; c2++) {
				uint8_t q[2] = {(uint8_t)c1, (uint8_t)c2};
				int64_t j;

				fmd_mem.n = bm_mem.n = 0;
				rb3_fmd_smem_TG(0, &fmi, 2, q, &fmd_mem, 1, 1);
				rb3_bmove_smem_TG(0, bm, 2, q, &bm_mem, 1, 1);
				if (fmd_mem.n != bm_mem.n) {
					fprintf(stderr, "FAIL: bmove_exhaustive 2-char [%d,%d]: n=%ld vs %ld\n",
						c1, c2, (long)fmd_mem.n, (long)bm_mem.n);
					ret = 1; goto done;
				}
				for (j = 0; j < (int64_t)fmd_mem.n; j++) {
					if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
						fprintf(stderr, "FAIL: bmove_exhaustive 2-char [%d,%d] mem[%ld] mismatch\n",
							c1, c2, (long)j);
						ret = 1; goto done;
					}
				}
				n_tested++;
			}
		}
	}

	/* All 3-char queries: 4^3 = 64 */
	{
		int c1, c2, c3;
		for (c1 = 1; c1 <= 4; c1++)
			for (c2 = 1; c2 <= 4; c2++)
				for (c3 = 1; c3 <= 4; c3++) {
					uint8_t q[3] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3};
					int64_t j;

					fmd_mem.n = bm_mem.n = 0;
					rb3_fmd_smem_TG(0, &fmi, 3, q, &fmd_mem, 1, 1);
					rb3_bmove_smem_TG(0, bm, 3, q, &bm_mem, 1, 1);
					if (fmd_mem.n != bm_mem.n) {
						fprintf(stderr, "FAIL: bmove_exhaustive 3-char [%d,%d,%d]: n=%ld vs %ld\n",
							c1, c2, c3, (long)fmd_mem.n, (long)bm_mem.n);
						ret = 1; goto done;
					}
					for (j = 0; j < (int64_t)fmd_mem.n; j++) {
						if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
							fprintf(stderr, "FAIL: bmove_exhaustive 3-char [%d,%d,%d] mem[%ld] mismatch\n",
								c1, c2, c3, (long)j);
							ret = 1; goto done;
						}
					}
					n_tested++;
				}
	}

	/* All 4-char queries: 4^4 = 256 */
	{
		int c1, c2, c3, c4;
		for (c1 = 1; c1 <= 4; c1++)
			for (c2 = 1; c2 <= 4; c2++)
				for (c3 = 1; c3 <= 4; c3++)
					for (c4 = 1; c4 <= 4; c4++) {
						uint8_t q[4] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3, (uint8_t)c4};
						int64_t j;

						/* Test both ORI and TG */
						fmd_mem.n = bm_mem.n = 0;
						rb3_fmd_smem(0, &fmi, 4, q, &fmd_mem, 1, 1);
						rb3_bmove_smem(0, bm, 4, q, &bm_mem, 1, 1);
						if (fmd_mem.n != bm_mem.n) {
							fprintf(stderr, "FAIL: bmove_exhaustive ORI 4-char: n mismatch\n");
							ret = 1; goto done;
						}
						for (j = 0; j < (int64_t)fmd_mem.n; j++) {
							if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
								fprintf(stderr, "FAIL: bmove_exhaustive ORI 4-char mem mismatch\n");
								ret = 1; goto done;
							}
						}

						fmd_mem.n = bm_mem.n = 0;
						rb3_fmd_smem_TG(0, &fmi, 4, q, &fmd_mem, 1, 1);
						rb3_bmove_smem_TG(0, bm, 4, q, &bm_mem, 1, 1);
						if (fmd_mem.n != bm_mem.n) {
							fprintf(stderr, "FAIL: bmove_exhaustive TG 4-char: n mismatch\n");
							ret = 1; goto done;
						}
						for (j = 0; j < (int64_t)fmd_mem.n; j++) {
							if (fmd_mem.a[j].x[0] != bm_mem.a[j].x[0] || fmd_mem.a[j].size != bm_mem.a[j].size || fmd_mem.a[j].info != bm_mem.a[j].info) {
								fprintf(stderr, "FAIL: bmove_exhaustive TG 4-char mem mismatch\n");
								ret = 1; goto done;
							}
						}
						n_tested++;
					}
	}

	fprintf(stderr, "test_bmove_smem_exhaustive: PASS (%d query patterns)\n", n_tested);
done:
	free(fmd_mem.a);
	free(bm_mem.a);
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test b-move init/destroy and edge cases.
 */
static int test_bmove_init(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_bmove_t *bm;
	int ret = 0;
	int c;
	int64_t i;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	bm = rb3_bmove_init(m);
	if (bm == 0) {
		fprintf(stderr, "FAIL: bmove_init returned NULL\n");
		ret = 1; goto done;
	}
	if (bm->mv != m) {
		fprintf(stderr, "FAIL: bmove_init mv pointer mismatch\n");
		ret = 1; goto done;
	}

	/* Verify cumrank: total should equal acc counts */
	for (c = 0; c < RB3_ASIZE; c++) {
		int64_t total = bm->cumrank[m->n_runs * RB3_ASIZE + c];
		int64_t expected = m->acc[c + 1] - m->acc[c];
		if (total != expected) {
			fprintf(stderr, "FAIL: bmove_init cumrank total[%d] = %ld, expected %ld\n",
				c, (long)total, (long)expected);
			ret = 1; goto done;
		}
	}

	/* Verify cumrank is monotonically non-decreasing for each character */
	for (c = 0; c < RB3_ASIZE; c++) {
		for (i = 1; i <= m->n_runs; i++) {
			if (bm->cumrank[i * RB3_ASIZE + c] < bm->cumrank[(i - 1) * RB3_ASIZE + c]) {
				fprintf(stderr, "FAIL: bmove_init cumrank[%ld][%d] not monotonic\n", (long)i, c);
				ret = 1; goto done;
			}
		}
	}

	/* Test that NULL input returns NULL */
	if (rb3_bmove_init(0) != 0) {
		fprintf(stderr, "FAIL: bmove_init(NULL) should return NULL\n");
		ret = 1; goto done;
	}

	/* Test destroy with NULL (should not crash) */
	rb3_bmove_destroy(0);

	fprintf(stderr, "test_bmove_init: PASS\n");
done:
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test that move_count produces the correct SA interval [lo, hi) at EACH step
 * of backward search, matching rank-based interval computation exactly.
 */
static int test_count_intervals(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	int ret = 0;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	rb3_move_precompute_dist(m);

	/* Build cumrank table (same as rb3_move_count does internally) */
	{
		int64_t *cumrank = (int64_t *)calloc(((int64_t)m->n_runs + 1) * RB3_ASIZE, sizeof(int64_t));
		int64_t j;
		int c2;
		for (j = 0; j < m->n_runs; j++) {
			for (c2 = 0; c2 < RB3_ASIZE; c2++)
				cumrank[(j + 1) * RB3_ASIZE + c2] = cumrank[j * RB3_ASIZE + c2];
			cumrank[(j + 1) * RB3_ASIZE + m->rows[j].c] += m->rows[j].len;
		}

		/* Test several multi-character patterns, verify intervals at each step */
		{
			int c1, c3;
			for (c1 = 0; c1 < RB3_ASIZE; c1++) {
				for (c2 = 0; c2 < RB3_ASIZE; c2++) {
					for (c3 = 0; c3 < RB3_ASIZE; c3++) {
						uint8_t pat[3] = {(uint8_t)c1, (uint8_t)c2, (uint8_t)c3};
						int64_t lo, hi, move_lo, move_hi;
						int step;

						/* Rank-based backward search, tracking interval at each step */
						lo = fmi.acc[pat[2]];
						hi = fmi.acc[pat[2] + 1];

						/* Move-based: same start */
						move_lo = lo;
						move_hi = hi;

						/* Verify initial interval matches */
						if (move_lo != lo || move_hi != hi) {
							fprintf(stderr, "FAIL: count_intervals [%d,%d,%d] step 0: move=[%ld,%ld), rank=[%ld,%ld)\n",
								c1, c2, c3, (long)move_lo, (long)move_hi, (long)lo, (long)hi);
							ret = 1; goto done;
						}

						for (step = 1; step >= 0 && lo < hi; --step) {
							int c = pat[step];
							int64_t ok[RB3_ASIZE], ol[RB3_ASIZE];

							/* Rank-based step */
							rb3_fmi_rank2a(&fmi, lo, hi, ok, ol);
							lo = fmi.acc[c] + ok[c];
							hi = fmi.acc[c] + ol[c];

							/* Move-based step using cumrank + binary search */
							{
								int64_t lo_run, hi_run, rank_lo, rank_hi;
								lo_run = 0;
								{
									int64_t a = 0, b = m->n_runs - 1;
									while (a < b) {
										int64_t mid = a + (b - a + 1) / 2;
										if (m->rows[mid].p <= move_lo) a = mid;
										else b = mid - 1;
									}
									lo_run = a;
								}
								if (move_hi < m->bwt_len) {
									int64_t a = 0, b = m->n_runs - 1;
									while (a < b) {
										int64_t mid = a + (b - a + 1) / 2;
										if (m->rows[mid].p <= move_hi) a = mid;
										else b = mid - 1;
									}
									hi_run = a;
								} else hi_run = m->n_runs - 1;

								rank_lo = cumrank[lo_run * RB3_ASIZE + c];
								if (m->rows[lo_run].c == c)
									rank_lo += move_lo - m->rows[lo_run].p;

								if (move_hi >= m->bwt_len)
									rank_hi = cumrank[m->n_runs * RB3_ASIZE + c];
								else {
									rank_hi = cumrank[hi_run * RB3_ASIZE + c];
									if (m->rows[hi_run].c == c)
										rank_hi += move_hi - m->rows[hi_run].p;
								}

								move_lo = m->acc[c] + rank_lo;
								move_hi = m->acc[c] + rank_hi;
							}

							if (move_lo != lo || move_hi != hi) {
								fprintf(stderr, "FAIL: count_intervals [%d,%d,%d] step %d: move=[%ld,%ld), rank=[%ld,%ld)\n",
									c1, c2, c3, 2 - step,
									(long)move_lo, (long)move_hi, (long)lo, (long)hi);
								free(cumrank);
								ret = 1; goto done;
							}
						}
					}
				}
			}
		}

		free(cumrank);
	}

	fprintf(stderr, "test_count_intervals: PASS (216 patterns, all steps verified)\n");
done:
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

/*
 * Test rank dispatch through fmi.bm: verify that when bm is set,
 * rb3_fmi_rank2a and rb3_fmi_rank1a produce the same results as direct FMD queries.
 */
static int test_rank_dispatch(void)
{
	uint8_t bwt[] = {2, 1, 1, 0, 2, 1, 4, 4, 1, 2};
	int64_t len = 10;
	rld_t *e;
	rb3_fmi_t fmi = {0};
	rb3_move_t *m;
	rb3_bmove_t *bm;
	int ret = 0;
	int64_t pos;

	e = rb3_enc_plain2rld(len, bwt, 3);
	rb3_fmi_init(&fmi, e, 0);
	m = rb3_move_build(&fmi);
	bm = rb3_bmove_init(m);

	/* First verify rank1a without bm (baseline) vs with bm for valid positions */
	for (pos = 0; pos < len; pos++) {
		int64_t ok_fmd[RB3_ASIZE], ok_bm[RB3_ASIZE];
		int c_fmd, c_bm;

		/* Direct FMD rank (bm not set) */
		c_fmd = rb3_fmi_rank1a(&fmi, pos, ok_fmd);

		/* Now set bm and use dispatch */
		fmi.bm = bm;
		c_bm = rb3_fmi_rank1a(&fmi, pos, ok_bm);
		fmi.bm = 0;

		if (c_fmd != c_bm) {
			fprintf(stderr, "FAIL: rank_dispatch rank1a pos=%ld: c_fmd=%d, c_bm=%d\n",
				(long)pos, c_fmd, c_bm);
			ret = 1; goto done;
		}
		{
			int c;
			for (c = 0; c < RB3_ASIZE; c++) {
				if (ok_fmd[c] != ok_bm[c]) {
					fprintf(stderr, "FAIL: rank_dispatch rank1a pos=%ld c=%d: fmd=%ld, bm=%ld\n",
						(long)pos, c, (long)ok_fmd[c], (long)ok_bm[c]);
					ret = 1; goto done;
				}
			}
		}
	}

	/* Verify rank2a: for all pairs [lo, hi) */
	{
		int64_t lo, hi;
		for (lo = 0; lo < len; lo++) {
			for (hi = lo + 1; hi <= len; hi++) {
				int64_t ok_fmd[RB3_ASIZE], ol_fmd[RB3_ASIZE];
				int64_t ok_bm[RB3_ASIZE], ol_bm[RB3_ASIZE];
				int c;

				rb3_fmi_rank2a(&fmi, lo, hi, ok_fmd, ol_fmd);
				fmi.bm = bm;
				rb3_fmi_rank2a(&fmi, lo, hi, ok_bm, ol_bm);
				fmi.bm = 0;

				for (c = 0; c < RB3_ASIZE; c++) {
					if (ok_fmd[c] != ok_bm[c] || ol_fmd[c] != ol_bm[c]) {
						fprintf(stderr, "FAIL: rank_dispatch rank2a [%ld,%ld) c=%d: fmd=(%ld,%ld), bm=(%ld,%ld)\n",
							(long)lo, (long)hi, c,
							(long)ok_fmd[c], (long)ol_fmd[c],
							(long)ok_bm[c], (long)ol_bm[c]);
						ret = 1; goto done;
					}
				}
			}
		}
	}

	/* Verify that backward search via rb3_fmi_extend1 works with bm dispatch */
	{
		int c;
		for (c = 0; c < RB3_ASIZE; c++) {
			int64_t k_fmd, l_fmd, k_bm, l_bm;

			k_fmd = fmi.acc[c];
			l_fmd = fmi.acc[c + 1];
			k_bm = k_fmd;
			l_bm = l_fmd;

			if (k_fmd < l_fmd) {
				int c2;
				for (c2 = 0; c2 < RB3_ASIZE; c2++) {
					int64_t sz_fmd, sz_bm;
					int64_t k1 = k_fmd, l1 = l_fmd, k2 = k_bm, l2 = l_bm;
					sz_fmd = rb3_fmi_extend1(&fmi, &k1, &l1, c2);
					fmi.bm = bm;
					sz_bm = rb3_fmi_extend1(&fmi, &k2, &l2, c2);
					fmi.bm = 0;
					if (sz_fmd != sz_bm || k1 != k2 || l1 != l2) {
						fprintf(stderr, "FAIL: rank_dispatch extend1 start=%d ext=%d: fmd=(%ld,%ld,%ld) bm=(%ld,%ld,%ld)\n",
							c, c2, (long)k1, (long)l1, (long)sz_fmd, (long)k2, (long)l2, (long)sz_bm);
						ret = 1; goto done;
					}
				}
			}
		}
	}

	fprintf(stderr, "test_rank_dispatch: PASS (rank1a, rank2a, extend1 all verified)\n");
done:
	fmi.bm = 0; /* detach before cleanup */
	rb3_bmove_destroy(bm);
	rb3_move_destroy(m);
	rb3_fmi_free(&fmi);
	return ret;
}

int main(void)
{
	int ret = 0;
	ret |= test_basic();
	ret |= test_single_char();
	ret |= test_alternating();
	ret |= test_fmr_backend();
	ret |= test_lf_all_positions();
	ret |= test_split_preserves_lf();
	ret |= test_split_d3();
	ret |= test_reposition();
	ret |= test_move_step();
	ret |= test_move_step_split();
	ret |= test_backward_search_chain();
	ret |= test_save_load();
	ret |= test_save_load_split();
	ret |= test_count_basic();
	ret |= test_count_larger();
	ret |= test_count_split();
	ret |= test_fmi_mv_field();
	ret |= test_bmove_init();
	ret |= test_bmove_extend();
	ret |= test_bmove_extend_split();
	ret |= test_bmove_smem();
	ret |= test_bmove_smem_params();
	ret |= test_bmove_smem_exhaustive();
	ret |= test_count_intervals();
	ret |= test_rank_dispatch();
	if (ret == 0)
		fprintf(stderr, "\nAll tests PASSED\n");
	else
		fprintf(stderr, "\nSome tests FAILED\n");
	return ret;
}
