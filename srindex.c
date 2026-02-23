#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "rb3priv.h"
#include "fm-index.h"
#include "srindex.h"
#include "rle.h"
#include "kalloc.h"
#include "kthread.h"
#include "ketopt.h"

/***************************
 * BWT run boundary scan   *
 ***************************/

/*
 * Scan the BWT to find all run boundaries. A run boundary is where the BWT
 * character changes. We record:
 *   - bwt_start[i]: the BWT position where run i starts
 *   - bwt_end[i]:   the BWT position of the last character in run i
 *
 * Returns the number of runs found.
 */
typedef struct {
	int64_t n, m;
	int64_t *bwt_start; /* BWT position of first char in each run */
	int64_t *bwt_end;   /* BWT position of last char in each run */
} run_bounds_t;

static run_bounds_t *scan_bwt_runs(const rb3_fmi_t *f)
{
	run_bounds_t *rb;
	int64_t pos = 0, l;
	int c, last_c = -1;

	rb = RB3_CALLOC(run_bounds_t, 1);
	rb->n = 0;
	rb->m = 1024;
	rb->bwt_start = RB3_MALLOC(int64_t, rb->m);
	rb->bwt_end = RB3_MALLOC(int64_t, rb->m);

	if (f->e) {
		rlditr_t itr;
		rld_itr_init(f->e, &itr, 0);
		while ((l = rld_dec(f->e, &itr, &c, 0)) > 0) {
			RB3_GROW(int64_t, rb->bwt_start, rb->n, rb->m);
			/* bwt_end must be grown in sync; re-grow if needed */
			if (rb->n >= (rb->m >> 1)) { /* ensure bwt_end has enough space too */
				int64_t old_m = rb->m;
				rb->m = rb->n + 1;
				rb->m += (rb->m >> 1) + 16;
				rb->bwt_start = RB3_REALLOC(int64_t, rb->bwt_start, rb->m);
				rb->bwt_end = RB3_REALLOC(int64_t, rb->bwt_end, rb->m);
				(void)old_m;
			}
			rb->bwt_start[rb->n] = pos;
			rb->bwt_end[rb->n] = pos + l - 1;
			rb->n++;
			pos += l;
		}
	} else if (f->r) {
		mritr_t ri;
		const uint8_t *block;
		mr_itr_first(f->r, &ri, 0);
		while ((block = mr_itr_next_block(&ri)) != 0) {
			const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
			while (q < end) {
				rle_dec1(q, c, l);
				if (c != last_c) {
					/* New run or continuation after block boundary */
					if (last_c >= 0 && rb->n > 0) {
						/* finalize the previous run's end if last_c changed */
					}
					RB3_GROW(int64_t, rb->bwt_start, rb->n, rb->m);
					if (rb->n >= (rb->m >> 1)) {
						rb->m = rb->n + 1;
						rb->m += (rb->m >> 1) + 16;
						rb->bwt_start = RB3_REALLOC(int64_t, rb->bwt_start, rb->m);
						rb->bwt_end = RB3_REALLOC(int64_t, rb->bwt_end, rb->m);
					}
					rb->bwt_start[rb->n] = pos;
					rb->bwt_end[rb->n] = pos + l - 1;
					rb->n++;
					last_c = c;
				} else {
					/* Same character, extend current run */
					rb->bwt_end[rb->n - 1] = pos + l - 1;
				}
				pos += l;
			}
		}
	}
	return rb;
}

static void run_bounds_destroy(run_bounds_t *rb)
{
	if (rb == 0) return;
	free(rb->bwt_start);
	free(rb->bwt_end);
	free(rb);
}

/***************************
 * SA computation at run   *
 * boundaries via backward *
 * walk from sentinels     *
 ***************************/

typedef struct {
	int64_t bwt_pos;
	int64_t sa_val;
} pos_sa_pair_t;

/* Sort helper for binary search on bwt_pos */
static int cmp_pos_sa(const void *a, const void *b)
{
	int64_t x = ((const pos_sa_pair_t*)a)->bwt_pos;
	int64_t y = ((const pos_sa_pair_t*)b)->bwt_pos;
	return (x > y) - (x < y);
}

/* Sort helper for int64_t */
static int cmp_int64(const void *a, const void *b)
{
	int64_t x = *(const int64_t*)a;
	int64_t y = *(const int64_t*)b;
	return (x > y) - (x < y);
}

/*
 * Walk backward from sentinel at BWT position k.
 *
 * Two-pass approach:
 * Pass 1: Walk to count total length.
 * Pass 2: Walk again, computing SA values on the fly. Record at target
 *          positions (run boundaries) and subsampled positions (SA % s == 0).
 */
static void sa_walk_one(const rb3_fmi_t *f, int64_t k, int32_t s,
                        const int64_t *targets, int64_t n_targets,
                        pos_sa_pair_t *tgt_out, int64_t *n_tgt_out,
                        pos_sa_pair_t *sub_out, int64_t max_sub, int64_t *n_sub_out,
                        int64_t *walk_dist_out, int64_t *dest_sent_out)
{
	int64_t ok[RB3_ASIZE];
	int32_t c;
	int64_t pos, dist, d, n_tgt = 0, n_sub = 0;

	/* Pass 1: count total walk length */
	pos = k;
	dist = 0;
	do {
		c = rb3_fmi_rank1a(f, pos, ok);
		pos = f->acc[c] + ok[c];
		dist++;
	} while (c);

	*walk_dist_out = dist;
	*dest_sent_out = pos; /* BWT position of destination sentinel (c was 0, pos = ok[0]) */

	/* Pass 2: walk again, recording */
	pos = k;
	d = 0;
	do {
		int64_t sa_val = dist - 1 - d;

		/* Check target (binary search) */
		{
			int64_t lo = 0, hi = n_targets;
			while (lo < hi) {
				int64_t mid = lo + (hi - lo) / 2;
				if (targets[mid] < pos) lo = mid + 1;
				else hi = mid;
			}
			if (lo < n_targets && targets[lo] == pos) {
				tgt_out[n_tgt].bwt_pos = pos;
				tgt_out[n_tgt].sa_val = sa_val;
				n_tgt++;
			}
		}

		/* Check subsampled (SA % s == 0) for s > 1 */
		if (s > 1 && sa_val % s == 0 && n_sub < max_sub) {
			sub_out[n_sub].bwt_pos = pos;
			sub_out[n_sub].sa_val = sa_val;
			n_sub++;
		}

		c = rb3_fmi_rank1a(f, pos, ok);
		pos = f->acc[c] + ok[c];
		d++;
	} while (c);

	*n_tgt_out = n_tgt;
	*n_sub_out = n_sub;
}

typedef struct {
	const rb3_fmi_t *f;
	const int64_t *targets;
	int64_t n_targets;
	int32_t s;
	/* Per-sentinel output buffers */
	pos_sa_pair_t **per_sent_tgt;
	int64_t *per_sent_tgt_n;
	pos_sa_pair_t **per_sent_sub;
	int64_t *per_sent_sub_n;
	int64_t total_n; /* total BWT length, for buffer sizing */
	int64_t n_sent;  /* number of sentinels */
	/* Per-sentinel walk info for multi-string correction */
	int64_t *walk_dist;  /* walk distance per sentinel */
	int64_t *dest_sent;  /* destination sentinel BWT position per sentinel */
} sa_mt_t;

static void sa_worker_func(void *data, long i, int tid)
{
	sa_mt_t *mt = (sa_mt_t*)data;
	int64_t max_sub;
	(void)tid;

	/* Estimate max subsampled entries per sentinel: n/(s*n_sent) + slack */
	if (mt->s > 1)
		max_sub = mt->total_n / mt->s + 2;
	else
		max_sub = 0;

	{
		/* Each sentinel walks ~n/n_sent positions; estimate targets proportionally.
		 * Use 4x the proportional estimate + slack, capped at n_targets. */
		int64_t max_tgt = mt->n_targets / (mt->n_sent > 0 ? mt->n_sent : 1) * 4 + 1024;
		if (max_tgt > mt->n_targets) max_tgt = mt->n_targets;
		mt->per_sent_tgt[i] = RB3_MALLOC(pos_sa_pair_t, max_tgt);
	}
	if (max_sub > 0)
		mt->per_sent_sub[i] = RB3_MALLOC(pos_sa_pair_t, max_sub);
	else
		mt->per_sent_sub[i] = 0;

	sa_walk_one(mt->f, (int64_t)i, mt->s,
	            mt->targets, mt->n_targets,
	            mt->per_sent_tgt[i], &mt->per_sent_tgt_n[i],
	            mt->per_sent_sub[i], max_sub, &mt->per_sent_sub_n[i],
	            &mt->walk_dist[i], &mt->dest_sent[i]);
}

/*
 * Compute SA values at all specified BWT positions and collect subsampled positions.
 */
static pos_sa_pair_t *compute_sa_at_positions(const rb3_fmi_t *f,
                                              int64_t *targets, int64_t n_targets,
                                              int32_t s, int n_threads,
                                              int64_t *n_results,
                                              pos_sa_pair_t **sub_results_out,
                                              int64_t *n_sub_out,
                                              int64_t **walk_dist_out,
                                              int64_t **dest_sent_out)
{
	sa_mt_t mt;
	int64_t n_sent = f->acc[1];
	int64_t i, total_tgt = 0, total_sub = 0;
	pos_sa_pair_t *all_tgt, *all_sub, *out;

	/* Sort targets for binary search */
	qsort(targets, n_targets, sizeof(int64_t), cmp_int64);

	mt.f = f;
	mt.targets = targets;
	mt.n_targets = n_targets;
	mt.s = s;
	mt.total_n = f->acc[RB3_ASIZE];
	mt.n_sent = n_sent;
	mt.per_sent_tgt = RB3_CALLOC(pos_sa_pair_t*, n_sent);
	mt.per_sent_tgt_n = RB3_CALLOC(int64_t, n_sent);
	mt.per_sent_sub = RB3_CALLOC(pos_sa_pair_t*, n_sent);
	mt.per_sent_sub_n = RB3_CALLOC(int64_t, n_sent);
	mt.walk_dist = RB3_CALLOC(int64_t, n_sent);
	mt.dest_sent = RB3_CALLOC(int64_t, n_sent);

	kt_for(n_threads, sa_worker_func, &mt, n_sent);

	/* Compute multi-string corrections.
	 * For multi-string BWTs, each sentinel walk produces per-sequence SA offsets.
	 * We need to convert these to absolute text positions by adding a correction
	 * equal to the cumulative length of preceding sequences.
	 * In a multi-string BWT, dest_sent[k] == k (self-loops) because each string
	 * is independent, so we use sequential sentinel order for corrections. */
	{
		int64_t *correction = RB3_CALLOC(int64_t, n_sent);
		if (n_sent > 1) {
			int64_t cum = 0;
			for (i = 0; i < n_sent; ++i) {
				correction[i] = cum;
				cum += mt.walk_dist[i];
			}
		}
		/* Apply corrections to target pairs */
		for (i = 0; i < n_sent; ++i) {
			int64_t j, corr = correction[i];
			for (j = 0; j < mt.per_sent_tgt_n[i]; ++j)
				mt.per_sent_tgt[i][j].sa_val += corr;
			for (j = 0; j < mt.per_sent_sub_n[i]; ++j)
				mt.per_sent_sub[i][j].sa_val += corr;
		}
		free(correction);
	}

	/* Merge target results */
	for (i = 0; i < n_sent; ++i)
		total_tgt += mt.per_sent_tgt_n[i];

	all_tgt = RB3_MALLOC(pos_sa_pair_t, total_tgt > 0 ? total_tgt : 1);
	out = all_tgt;
	for (i = 0; i < n_sent; ++i) {
		if (mt.per_sent_tgt_n[i] > 0)
			memcpy(out, mt.per_sent_tgt[i], mt.per_sent_tgt_n[i] * sizeof(pos_sa_pair_t));
		out += mt.per_sent_tgt_n[i];
		free(mt.per_sent_tgt[i]);
	}
	qsort(all_tgt, total_tgt, sizeof(pos_sa_pair_t), cmp_pos_sa);

	/* Merge subsampled results */
	if (s > 1) {
		for (i = 0; i < n_sent; ++i)
			total_sub += mt.per_sent_sub_n[i];
		all_sub = RB3_MALLOC(pos_sa_pair_t, total_sub > 0 ? total_sub : 1);
		out = all_sub;
		for (i = 0; i < n_sent; ++i) {
			if (mt.per_sent_sub_n[i] > 0)
				memcpy(out, mt.per_sent_sub[i], mt.per_sent_sub_n[i] * sizeof(pos_sa_pair_t));
			out += mt.per_sent_sub_n[i];
			free(mt.per_sent_sub[i]);
		}
		qsort(all_sub, total_sub, sizeof(pos_sa_pair_t), cmp_pos_sa);
	} else {
		all_sub = 0;
		for (i = 0; i < n_sent; ++i)
			free(mt.per_sent_sub[i]);
	}

	free(mt.per_sent_tgt);
	free(mt.per_sent_tgt_n);
	free(mt.per_sent_sub);
	free(mt.per_sent_sub_n);

	*n_results = total_tgt;
	*sub_results_out = all_sub;
	*n_sub_out = total_sub;
	*walk_dist_out = mt.walk_dist;
	*dest_sent_out = mt.dest_sent;
	return all_tgt;
}

/***************************
 * Phi function building   *
 ***************************/

/* Sort helper: sort by SA value */
static int cmp_pos_sa_by_sa(const void *a, const void *b)
{
	int64_t x = ((const pos_sa_pair_t*)a)->sa_val;
	int64_t y = ((const pos_sa_pair_t*)b)->sa_val;
	return (x > y) - (x < y);
}

static void build_phi(rb3_srindex_t *sr, const run_bounds_t *rb,
                      const pos_sa_pair_t *sa_pairs, int64_t n_pairs)
{
	int64_t i;
	pos_sa_pair_t *start_pairs;

	sr->n_runs = rb->n;
	sr->phi_sa = RB3_MALLOC(int64_t, rb->n);
	sr->phi_da = RB3_MALLOC(int64_t, rb->n);

	/* For each run start, look up SA value */
	start_pairs = RB3_MALLOC(pos_sa_pair_t, rb->n);
	for (i = 0; i < rb->n; ++i) {
		int64_t lo = 0, hi = n_pairs, target = rb->bwt_start[i];
		while (lo < hi) {
			int64_t mid = lo + (hi - lo) / 2;
			if (sa_pairs[mid].bwt_pos < target) lo = mid + 1;
			else hi = mid;
		}
		assert(lo < n_pairs && sa_pairs[lo].bwt_pos == target);
		start_pairs[i].bwt_pos = target;
		start_pairs[i].sa_val = sa_pairs[lo].sa_val;
	}

	/* Collect (SA_at_start, SA_at_prev_position) pairs */
	for (i = 0; i < rb->n; ++i) {
		int64_t prev_pos;
		sr->phi_sa[i] = start_pairs[i].sa_val;

		if (rb->bwt_start[i] == 0) {
			sr->phi_da[i] = -1;
			continue;
		}

		prev_pos = rb->bwt_start[i] - 1;
		{
			int64_t lo = 0, hi = n_pairs;
			while (lo < hi) {
				int64_t mid = lo + (hi - lo) / 2;
				if (sa_pairs[mid].bwt_pos < prev_pos) lo = mid + 1;
				else hi = mid;
			}
			assert(lo < n_pairs && sa_pairs[lo].bwt_pos == prev_pos);
			sr->phi_da[i] = sa_pairs[lo].sa_val;
		}
	}

	free(start_pairs);

	/* Sort phi_sa[] and phi_da[] together by phi_sa value for binary search */
	{
		pos_sa_pair_t *tmp = RB3_MALLOC(pos_sa_pair_t, rb->n);
		for (i = 0; i < rb->n; ++i) {
			tmp[i].bwt_pos = sr->phi_da[i]; /* reuse bwt_pos field for phi_da */
			tmp[i].sa_val = sr->phi_sa[i];
		}
		qsort(tmp, rb->n, sizeof(pos_sa_pair_t), cmp_pos_sa_by_sa);
		for (i = 0; i < rb->n; ++i) {
			sr->phi_sa[i] = tmp[i].sa_val;
			sr->phi_da[i] = tmp[i].bwt_pos;
		}
		free(tmp);
	}
}

/***************************
 * Toehold support         *
 ***************************/

static void build_toehold(rb3_srindex_t *sr, const run_bounds_t *rb,
                          const pos_sa_pair_t *sa_pairs, int64_t n_pairs)
{
	int64_t i;

	sr->run_pos = RB3_MALLOC(int64_t, rb->n);
	sr->run_sa = RB3_MALLOC(int64_t, rb->n);
	sr->n_samples = rb->n;

	for (i = 0; i < rb->n; ++i) {
		int64_t target = rb->bwt_end[i];
		int64_t lo = 0, hi = n_pairs;
		while (lo < hi) {
			int64_t mid = lo + (hi - lo) / 2;
			if (sa_pairs[mid].bwt_pos < target) lo = mid + 1;
			else hi = mid;
		}
		assert(lo < n_pairs && sa_pairs[lo].bwt_pos == target);
		sr->run_pos[i] = target;
		sr->run_sa[i] = sa_pairs[lo].sa_val;
	}
}

/***************************
 * Subsampled SA building  *
 ***************************/

/* Build bitvector marking sub_pos positions for O(1) lookup in locate_one */
static void build_sub_bitvector(rb3_srindex_t *sr)
{
	int64_t i, n_words = (sr->n + 63) / 64;
	sr->sub_bv = RB3_CALLOC(uint64_t, n_words > 0 ? n_words : 1);
	for (i = 0; i < sr->n_sub; ++i) {
		int64_t p = sr->sub_pos[i];
		if (p >= 0 && p < sr->n)
			sr->sub_bv[p >> 6] |= 1ULL << (p & 63);
	}
}

static void build_subsampled(rb3_srindex_t *sr, int32_t s,
                             const pos_sa_pair_t *sub_pairs, int64_t n_sub)
{
	int64_t i;
	if (s <= 1) {
		/* For s=1, alias run boundary samples (no copy needed) */
		sr->n_sub = sr->n_samples;
		sr->sub_pos = sr->run_pos;
		sr->sub_sa = sr->run_sa;
		sr->sub_is_alias = 1;
	} else {
		sr->n_sub = n_sub;
		sr->sub_pos = RB3_MALLOC(int64_t, n_sub > 0 ? n_sub : 1);
		sr->sub_sa = RB3_MALLOC(int64_t, n_sub > 0 ? n_sub : 1);
		for (i = 0; i < n_sub; ++i) {
			sr->sub_pos[i] = sub_pairs[i].bwt_pos;
			sr->sub_sa[i] = sub_pairs[i].sa_val;
		}
		sr->sub_is_alias = 0;
	}
	build_sub_bitvector(sr);
}

/***************************
 * Public API              *
 ***************************/

rb3_srindex_t *rb3_srindex_build(const void *f_, int32_t s, int n_threads)
{
	const rb3_fmi_t *f = (const rb3_fmi_t*)f_;
	rb3_srindex_t *sr;
	run_bounds_t *rb;
	int64_t *targets, n_targets, n_pairs, n_sub, i;
	int64_t *walk_dist, *dest_sent;
	pos_sa_pair_t *sa_pairs, *sub_pairs;

	if (s < 1) s = 1;
	if (f->e == 0 && f->r == 0) return 0;

	/* Step 1: Scan BWT to find run boundaries */
	rb = scan_bwt_runs(f);
	if (rb == 0 || rb->n == 0) {
		run_bounds_destroy(rb);
		return 0;
	}

	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] found %lld BWT runs, s=%d\n", __func__, (long long)rb->n, s);

	/* Step 2: Collect target BWT positions (run starts + ends) */
	n_targets = 0;
	targets = RB3_MALLOC(int64_t, rb->n * 2);
	for (i = 0; i < rb->n; ++i) {
		targets[n_targets++] = rb->bwt_start[i];
		targets[n_targets++] = rb->bwt_end[i];
	}
	/* Deduplicate */
	qsort(targets, n_targets, sizeof(int64_t), cmp_int64);
	{
		int64_t j = 0;
		for (i = 0; i < n_targets; ++i)
			if (i == 0 || targets[i] != targets[i-1])
				targets[j++] = targets[i];
		n_targets = j;
	}

	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] computing SA at %lld BWT positions\n", __func__, (long long)n_targets);

	/* Step 3: Compute SA values at targets + collect subsampled positions */
	sa_pairs = compute_sa_at_positions(f, targets, n_targets, s, n_threads,
	                                   &n_pairs, &sub_pairs, &n_sub,
	                                   &walk_dist, &dest_sent);
	free(targets);

	if (rb3_verbose >= 3) {
		fprintf(stderr, "[M::%s] computed %lld SA values", __func__, (long long)n_pairs);
		if (s > 1)
			fprintf(stderr, ", %lld subsampled positions (s=%d)", (long long)n_sub, s);
		fprintf(stderr, "\n");
	}

	/* Step 4: Build the SR-index */
	sr = RB3_CALLOC(rb3_srindex_t, 1);
	sr->n = f->acc[RB3_ASIZE];
	sr->s = s;
	sr->m = f->acc[1];

	build_phi(sr, rb, sa_pairs, n_pairs);
	build_toehold(sr, rb, sa_pairs, n_pairs);
	build_subsampled(sr, s, sub_pairs, n_sub);

	/* Step 5: Build multi-string mapping (cum_len and text_order_sid).
	 * In a multi-string BWT, dest_sent[k] == k (self-loops) because each
	 * string is independent. We use sequential sentinel order: sentinel k
	 * maps to text_order_sid[k] = k, and cum_len[k] = sum of walk_dist[0..k-1]. */
	{
		int64_t m = sr->m, cum = 0;
		sr->cum_len = RB3_MALLOC(int64_t, m + 1);
		sr->text_order_sid = RB3_MALLOC(int64_t, m > 0 ? m : 1);
		for (i = 0; i < m; ++i) {
			sr->text_order_sid[i] = i;
			sr->cum_len[i] = cum;
			cum += walk_dist[i];
		}
		sr->cum_len[m] = cum;
	}

	free(sa_pairs);
	free(sub_pairs);
	free(walk_dist);
	free(dest_sent);
	run_bounds_destroy(rb);

	return sr;
}

int64_t rb3_srindex_phi(const rb3_srindex_t *sr, int64_t sa_val)
{
	int64_t lo = 0, hi = sr->n_runs;

	if (sr == 0 || sr->n_runs == 0) return -1;

	/* Binary search for the largest i such that phi_sa[i] <= sa_val */
	while (lo < hi) {
		int64_t mid = lo + (hi - lo) / 2;
		if (sr->phi_sa[mid] <= sa_val) lo = mid + 1;
		else hi = mid;
	}
	if (lo == 0) return -1;
	lo--;

	return sr->phi_da[lo] + (sa_val - sr->phi_sa[lo]);
}

int64_t rb3_srindex_toehold(const rb3_srindex_t *sr, int64_t bwt_pos)
{
	int64_t lo = 0, hi;

	if (sr == 0 || sr->n_samples == 0) return -1;
	hi = sr->n_samples;

	while (lo < hi) {
		int64_t mid = lo + (hi - lo) / 2;
		if (sr->run_pos[mid] < bwt_pos) lo = mid + 1;
		else hi = mid;
	}
	if (lo < sr->n_samples && sr->run_pos[lo] == bwt_pos)
		return sr->run_sa[lo];
	return -1;
}

int64_t rb3_srindex_locate(const rb3_srindex_t *sr, int64_t lo, int64_t hi,
                           int64_t toehold_sa, int64_t *out)
{
	int64_t i, n = hi - lo;

	if (n <= 0) return 0;

	out[n - 1] = toehold_sa;
	for (i = n - 2; i >= 0; --i) {
		out[i] = rb3_srindex_phi(sr, out[i + 1]);
		if (out[i] < 0) return -1;
	}
	return n;
}

int64_t rb3_srindex_locate_one(const rb3_srindex_t *sr, const void *f_, int64_t bwt_pos)
{
	const rb3_fmi_t *f = (const rb3_fmi_t*)f_;
	int64_t ok[RB3_ASIZE];
	int32_t c;
	int64_t pos = bwt_pos, steps = 0;

	if (sr == 0 || sr->n_sub == 0) return -1;

	/* Walk LF from bwt_pos, checking subsampled set at each step.
	 * SA[LF(pos)] = SA[pos] - 1, so after j steps SA = SA[bwt_pos] - j.
	 * We stop when we find a stored sample. Return: stored_sa + steps. */
	while (steps <= sr->s + sr->n) { /* safety bound */
		/* O(1) bitvector test + O(log n_sub) binary search only on hit */
		if (sr->sub_bv && pos >= 0 && pos < sr->n &&
		    (sr->sub_bv[pos >> 6] & (1ULL << (pos & 63)))) {
			/* Confirmed in bitvector; binary search for exact SA value */
			int64_t lo = 0, hi = sr->n_sub;
			while (lo < hi) {
				int64_t mid = lo + (hi - lo) / 2;
				if (sr->sub_pos[mid] < pos) lo = mid + 1;
				else hi = mid;
			}
			if (lo < sr->n_sub && sr->sub_pos[lo] == pos)
				return sr->sub_sa[lo] + steps;
		}

		/* LF step */
		c = rb3_fmi_rank1a(f, pos, ok);
		pos = f->acc[c] + ok[c];
		steps++;
		if (c == 0) {
			if (sr->cum_len && pos < sr->m)
				return sr->cum_len[pos] + (steps - 1);
			break;
		}
	}
	return -1;
}

int64_t rb3_srindex_locate_all(const rb3_srindex_t *sr, const void *f_,
                               int64_t lo, int64_t hi, int64_t *positions, int64_t max_pos)
{
	int64_t n, i, toehold_sa;

	if (sr == 0) return -1;
	n = hi - lo;
	if (n <= 0) return 0;
	if (n > max_pos) n = max_pos;

	/* Resolve toehold SA[hi-1]:
	 * First try direct toehold lookup (works if hi-1 is at a run boundary) */
	toehold_sa = rb3_srindex_toehold(sr, hi - 1);
	if (toehold_sa < 0) {
		/* hi-1 is not at a run boundary; use LF walking */
		toehold_sa = rb3_srindex_locate_one(sr, f_, hi - 1);
	}
	if (toehold_sa < 0) return -1;

	/* Enumerate using phi: SA[hi-1], SA[hi-2], ..., SA[lo] */
	positions[n - 1] = toehold_sa;
	for (i = n - 2; i >= 0; --i) {
		positions[i] = rb3_srindex_phi(sr, positions[i + 1]);
		if (positions[i] < 0) return -1;
	}
	return n;
}

int64_t rb3_srindex_multi(void *km, const rb3_fmi_t *f, const rb3_srindex_t *sr,
                          int64_t lo, int64_t hi, int64_t max_pos, rb3_pos_t *pos)
{
	int64_t i, n, *sa_vals;
	n = hi - lo;
	if (n <= 0) return 0;
	if (n > max_pos) n = max_pos;
	sa_vals = (int64_t*)kmalloc(km, n * sizeof(int64_t));
	n = rb3_srindex_locate_all(sr, f, lo, hi, sa_vals, n);
	if (n < 0) { kfree(km, sa_vals); return 0; }
	for (i = 0; i < n; ++i) {
		int64_t sa = sa_vals[i];
		int64_t lo2 = 0, hi2 = sr->m;
		while (lo2 < hi2) {
			int64_t mid = lo2 + (hi2 - lo2) / 2;
			if (sr->cum_len[mid + 1] <= sa) lo2 = mid + 1;
			else hi2 = mid;
		}
		pos[i].sid = sr->text_order_sid[lo2];
		pos[i].pos = sa - sr->cum_len[lo2];
	}
	kfree(km, sa_vals);
	return n;
}

void rb3_srindex_destroy(rb3_srindex_t *sr)
{
	if (sr == 0) return;
	free(sr->phi_sa);
	free(sr->phi_da);
	free(sr->run_pos);
	free(sr->run_sa);
	if (!sr->sub_is_alias) {
		free(sr->sub_pos);
		free(sr->sub_sa);
	}
	free(sr->sub_bv);
	free(sr->cum_len);
	free(sr->text_order_sid);
	free(sr);
}

/***************************
 * SR-index serialization  *
 ***************************/

/*
 * V3 format: compressed SR-index with three optimizations:
 *
 * 1. 32-bit integer mode: For n < 2^32, sorted array samples and all
 *    absolute values stored as uint32 instead of int64. Halves the index.
 *
 * 2. Delta-encoded sorted arrays (phi_sa, run_pos, sub_pos):
 *    - Absolute sample every K entries (K = DELTA_SAMPLE_K)
 *    - Between samples: 16-bit or 32-bit deltas (chosen per array)
 *    - Provides O(log(count/K) + K) binary search via sample index
 *
 * 3. Bit-packed unsorted arrays (phi_da, run_sa, sub_sa):
 *    - ceil(log2(n)) bits per entry packed into a byte array
 *    - Access via (index * bits) >> 3 shift+mask
 *
 * Header (52 bytes):
 *   magic "SRI\3"     4 bytes
 *   s                  4 bytes (int32)
 *   m                  8 bytes (int64)
 *   n                  8 bytes (int64)
 *   n_runs             8 bytes (int64)
 *   n_samples          8 bytes (int64)
 *   n_sub              8 bytes (int64) - 0 if s<=1 (alias)
 *   bit_width          1 byte  - ceil(log2(n)) for bit-packing
 *   delta_bits         1 byte  - 16 or 32 for delta encoding
 *   reserved           2 bytes
 *
 * Array layout:
 *   phi_sa     delta-encoded
 *   phi_da     bit-packed
 *   run_pos    delta-encoded
 *   run_sa     bit-packed
 *   sub_pos    delta-encoded (if not alias)
 *   sub_sa     bit-packed (if not alias)
 *   cum_len    raw int64 (small)
 *   tosid      raw int64 (small)
 */

#define DELTA_SAMPLE_K 64

/* Compute ceil(log2(n)), minimum 1 */
static int compute_bit_width(int64_t n)
{
	int bits = 0;
	int64_t v = n - 1;
	if (v <= 0) return 1;
	while (v > 0) { bits++; v >>= 1; }
	return bits;
}

/* Check if any delta in a sorted array exceeds 16-bit range */
static int need_32bit_deltas(const int64_t *arr, int64_t count)
{
	int64_t i;
	for (i = 1; i < count; ++i)
		if (arr[i] - arr[i-1] > 65535) return 1;
	return 0;
}

/*
 * Write a delta-encoded sorted array.
 * Format: [n_samples uint32 absolute samples] [count deltas of delta_bits/8 bytes each]
 * Sample at positions 0, K, 2K, ...
 */
static void write_delta_array(FILE *fp, const int64_t *arr, int64_t count,
                               int delta_bits)
{
	int64_t i, n_samp;
	if (count == 0) return;

	/* Write absolute samples every K entries */
	n_samp = (count + DELTA_SAMPLE_K - 1) / DELTA_SAMPLE_K;
	for (i = 0; i < n_samp; ++i) {
		uint32_t v = (uint32_t)arr[i * DELTA_SAMPLE_K];
		fwrite(&v, 4, 1, fp);
	}

	/* Write deltas for all entries after the first */
	if (delta_bits == 16) {
		for (i = 1; i < count; ++i) {
			uint16_t d = (uint16_t)(arr[i] - arr[i-1]);
			fwrite(&d, 2, 1, fp);
		}
	} else {
		for (i = 1; i < count; ++i) {
			uint32_t d = (uint32_t)(arr[i] - arr[i-1]);
			fwrite(&d, 4, 1, fp);
		}
	}
}

/*
 * Read a delta-encoded sorted array back into int64_t.
 */
static void read_delta_array(FILE *fp, int64_t *arr, int64_t count,
                              int delta_bits)
{
	int64_t i, n_samp;
	uint32_t *samples;
	if (count == 0) return;

	n_samp = (count + DELTA_SAMPLE_K - 1) / DELTA_SAMPLE_K;
	samples = RB3_MALLOC(uint32_t, n_samp);
	fread(samples, 4, n_samp, fp);

	/* Set sample positions */
	for (i = 0; i < n_samp; ++i)
		arr[i * DELTA_SAMPLE_K] = samples[i];
	free(samples);

	/* Read and accumulate deltas */
	if (delta_bits == 16) {
		for (i = 1; i < count; ++i) {
			uint16_t d;
			fread(&d, 2, 1, fp);
			if (i % DELTA_SAMPLE_K != 0) /* don't overwrite samples */
				arr[i] = arr[i-1] + d;
		}
	} else {
		for (i = 1; i < count; ++i) {
			uint32_t d;
			fread(&d, 4, 1, fp);
			if (i % DELTA_SAMPLE_K != 0)
				arr[i] = arr[i-1] + d;
		}
	}
}

/*
 * Write a bit-packed array. Each value uses 'bits' bits, packed into bytes.
 * Total bytes: ceil(count * bits / 8)
 */
static void write_packed_array(FILE *fp, const int64_t *arr, int64_t count,
                                int bits)
{
	int64_t i, n_bytes;
	uint8_t *buf;
	if (count == 0 || bits == 0) return;

	n_bytes = ((int64_t)count * bits + 7) / 8;
	buf = RB3_CALLOC(uint8_t, n_bytes);

	for (i = 0; i < count; ++i) {
		uint64_t val = (uint64_t)arr[i];
		int64_t bit_pos = i * bits;
		int64_t byte_pos = bit_pos >> 3;
		int bit_off = (int)(bit_pos & 7);
		int remaining = bits;

		/* Write value across byte boundaries */
		while (remaining > 0) {
			int can_write = 8 - bit_off;
			if (can_write > remaining) can_write = remaining;
			buf[byte_pos] |= (uint8_t)((val & ((1U << can_write) - 1)) << bit_off);
			val >>= can_write;
			remaining -= can_write;
			byte_pos++;
			bit_off = 0;
		}
	}
	fwrite(buf, 1, n_bytes, fp);
	free(buf);
}

/*
 * Read a bit-packed array back into int64_t.
 */
static void read_packed_array(FILE *fp, int64_t *arr, int64_t count,
                               int bits)
{
	int64_t i, n_bytes;
	uint8_t *buf;
	if (count == 0 || bits == 0) return;

	n_bytes = ((int64_t)count * bits + 7) / 8;
	buf = RB3_MALLOC(uint8_t, n_bytes);
	fread(buf, 1, n_bytes, fp);

	for (i = 0; i < count; ++i) {
		uint64_t val = 0;
		int64_t bit_pos = i * bits;
		int64_t byte_pos = bit_pos >> 3;
		int bit_off = (int)(bit_pos & 7);
		int remaining = bits;
		int shift = 0;

		while (remaining > 0) {
			int can_read = 8 - bit_off;
			if (can_read > remaining) can_read = remaining;
			val |= (uint64_t)((buf[byte_pos] >> bit_off) & ((1U << can_read) - 1)) << shift;
			shift += can_read;
			remaining -= can_read;
			byte_pos++;
			bit_off = 0;
		}
		arr[i] = (int64_t)val;
	}
	free(buf);
}

int rb3_srindex_dump(const rb3_srindex_t *sr, const char *fn)
{
	FILE *fp;
	int32_t y;
	int64_t n_sub_disk;
	int bit_width, delta_bits;
	uint8_t hdr_extra[4];

	if (sr == 0) return -1;
	fp = fn && strcmp(fn, "-") ? fopen(fn, "wb") : fdopen(1, "wb");
	if (fp == 0) return -1;

	bit_width = compute_bit_width(sr->n);
	delta_bits = (need_32bit_deltas(sr->phi_sa, sr->n_runs) ||
	              need_32bit_deltas(sr->run_pos, sr->n_samples) ||
	              (!sr->sub_is_alias && need_32bit_deltas(sr->sub_pos, sr->n_sub)))
	             ? 32 : 16;

	/*
	 * Format version 3: compressed with delta-encoding and bit-packing.
	 * For s<=1, n_sub is written as 0 to indicate alias.
	 */
	fwrite("SRI\3", 1, 4, fp);
	y = sr->s; fwrite(&y, 4, 1, fp);
	fwrite(&sr->m, 8, 1, fp);
	fwrite(&sr->n, 8, 1, fp);
	fwrite(&sr->n_runs, 8, 1, fp);
	fwrite(&sr->n_samples, 8, 1, fp);
	n_sub_disk = sr->sub_is_alias ? 0 : sr->n_sub;
	fwrite(&n_sub_disk, 8, 1, fp);
	hdr_extra[0] = (uint8_t)bit_width;
	hdr_extra[1] = (uint8_t)delta_bits;
	hdr_extra[2] = 0;
	hdr_extra[3] = 0;
	fwrite(hdr_extra, 1, 4, fp);

	/* Sorted arrays: delta-encoded */
	write_delta_array(fp, sr->phi_sa, sr->n_runs, delta_bits);
	/* Unsorted: phi_da can contain -1, so we map -1 -> (1<<bw)-1 before packing.
	 * Use bw+1 bits to avoid collision when n-1 == (1<<bw)-1. */
	{
		int64_t ii;
		int phi_da_bits = bit_width + 1;
		int64_t *phi_da_mapped = RB3_MALLOC(int64_t, sr->n_runs > 0 ? sr->n_runs : 1);
		for (ii = 0; ii < sr->n_runs; ++ii)
			phi_da_mapped[ii] = sr->phi_da[ii] < 0 ? ((1LL << phi_da_bits) - 1) : sr->phi_da[ii];
		write_packed_array(fp, phi_da_mapped, sr->n_runs, phi_da_bits);
		free(phi_da_mapped);
	}
	/* Sorted arrays: delta-encoded */
	write_delta_array(fp, sr->run_pos, sr->n_samples, delta_bits);
	/* Unsorted arrays: bit-packed */
	write_packed_array(fp, sr->run_sa, sr->n_samples, bit_width);

	if (!sr->sub_is_alias) {
		write_delta_array(fp, sr->sub_pos, sr->n_sub, delta_bits);
		write_packed_array(fp, sr->sub_sa, sr->n_sub, bit_width);
	}

	/* Small arrays: raw int64 */
	fwrite(sr->cum_len, 8, sr->m + 1, fp);
	fwrite(sr->text_order_sid, 8, sr->m, fp);
	fclose(fp);
	return 0;
}

rb3_srindex_t *rb3_srindex_restore(const char *fn)
{
	FILE *fp;
	int32_t y, version;
	char magic[4];
	rb3_srindex_t *sr;

	fp = fn && strcmp(fn, "-") ? fopen(fn, "rb") : fdopen(0, "rb");
	if (fp == 0) return 0;
	if (fread(magic, 1, 4, fp) != 4 || magic[0] != 'S' || magic[1] != 'R' || magic[2] != 'I') {
		fclose(fp);
		return 0;
	}
	version = (unsigned char)magic[3];
	if (version < 1 || version > 3) {
		fclose(fp);
		return 0;
	}
	sr = RB3_CALLOC(rb3_srindex_t, 1);
	fread(&y, 4, 1, fp); sr->s = y;
	fread(&sr->m, 8, 1, fp);
	fread(&sr->n, 8, 1, fp);
	fread(&sr->n_runs, 8, 1, fp);
	fread(&sr->n_samples, 8, 1, fp);
	fread(&sr->n_sub, 8, 1, fp);

	if (version == 3) {
		uint8_t hdr_extra[4];
		int bit_width, delta_bits;

		fread(hdr_extra, 1, 4, fp);
		bit_width = hdr_extra[0];
		delta_bits = hdr_extra[1];

		sr->phi_sa = RB3_MALLOC(int64_t, sr->n_runs > 0 ? sr->n_runs : 1);
		sr->phi_da = RB3_MALLOC(int64_t, sr->n_runs > 0 ? sr->n_runs : 1);
		sr->run_pos = RB3_MALLOC(int64_t, sr->n_samples > 0 ? sr->n_samples : 1);
		sr->run_sa = RB3_MALLOC(int64_t, sr->n_samples > 0 ? sr->n_samples : 1);

		read_delta_array(fp, sr->phi_sa, sr->n_runs, delta_bits);
		/* phi_da uses bit_width+1 bits; sentinel (all 1s) maps back to -1 */
		{
			int phi_da_bits = bit_width + 1;
			int64_t sentinel = (1LL << phi_da_bits) - 1;
			int64_t ii;
			read_packed_array(fp, sr->phi_da, sr->n_runs, phi_da_bits);
			for (ii = 0; ii < sr->n_runs; ++ii)
				if (sr->phi_da[ii] == sentinel) sr->phi_da[ii] = -1;
		}
		read_delta_array(fp, sr->run_pos, sr->n_samples, delta_bits);
		read_packed_array(fp, sr->run_sa, sr->n_samples, bit_width);

		if (sr->n_sub == 0 && sr->s <= 1) {
			sr->n_sub = sr->n_samples;
			sr->sub_pos = sr->run_pos;
			sr->sub_sa = sr->run_sa;
			sr->sub_is_alias = 1;
		} else {
			sr->sub_pos = RB3_MALLOC(int64_t, sr->n_sub > 0 ? sr->n_sub : 1);
			sr->sub_sa = RB3_MALLOC(int64_t, sr->n_sub > 0 ? sr->n_sub : 1);
			read_delta_array(fp, sr->sub_pos, sr->n_sub, delta_bits);
			read_packed_array(fp, sr->sub_sa, sr->n_sub, bit_width);
			sr->sub_is_alias = 0;
		}
	} else {
		/* v1/v2: raw int64 arrays */
		sr->phi_sa = RB3_MALLOC(int64_t, sr->n_runs > 0 ? sr->n_runs : 1);
		sr->phi_da = RB3_MALLOC(int64_t, sr->n_runs > 0 ? sr->n_runs : 1);
		sr->run_pos = RB3_MALLOC(int64_t, sr->n_samples > 0 ? sr->n_samples : 1);
		sr->run_sa = RB3_MALLOC(int64_t, sr->n_samples > 0 ? sr->n_samples : 1);
		fread(sr->phi_sa, 8, sr->n_runs, fp);
		fread(sr->phi_da, 8, sr->n_runs, fp);
		fread(sr->run_pos, 8, sr->n_samples, fp);
		fread(sr->run_sa, 8, sr->n_samples, fp);

		if (version >= 2 && sr->n_sub == 0 && sr->s <= 1) {
			sr->n_sub = sr->n_samples;
			sr->sub_pos = sr->run_pos;
			sr->sub_sa = sr->run_sa;
			sr->sub_is_alias = 1;
		} else {
			sr->sub_pos = RB3_MALLOC(int64_t, sr->n_sub > 0 ? sr->n_sub : 1);
			sr->sub_sa = RB3_MALLOC(int64_t, sr->n_sub > 0 ? sr->n_sub : 1);
			fread(sr->sub_pos, 8, sr->n_sub, fp);
			fread(sr->sub_sa, 8, sr->n_sub, fp);
			sr->sub_is_alias = 0;
		}
	}

	sr->cum_len = RB3_MALLOC(int64_t, sr->m + 1);
	sr->text_order_sid = RB3_MALLOC(int64_t, sr->m > 0 ? sr->m : 1);
	fread(sr->cum_len, 8, sr->m + 1, fp);
	fread(sr->text_order_sid, 8, sr->m, fp);
	fclose(fp);

	/* Build bitvector for fast locate_one */
	build_sub_bitvector(sr);

	return sr;
}

/***************************
 * main_srindex CLI        *
 ***************************/

int main_srindex(int argc, char *argv[])
{
	int c, n_threads = 4, s_param = 8;
	rb3_srindex_t *sr;
	rb3_fmi_t f;
	char *fn = 0;
	ketopt_t o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, "t:s:o:", 0)) >= 0) {
		if (c == 't') n_threads = atoi(o.arg);
		else if (c == 's') s_param = atoi(o.arg);
		else if (c == 'o') fn = o.arg;
	}
	if (argc == o.ind) {
		fprintf(stderr, "Usage: ropebwt3 srindex [options] <in.fmd>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "  -s INT     subsampling parameter [%d]\n", s_param);
		fprintf(stderr, "  -o FILE    output file [<in.fmd>.sri]\n");
		return 1;
	}
	rb3_fmi_restore(&f, argv[o.ind], 0);
	if (f.e == 0 && f.r == 0) {
		fprintf(stderr, "[E::%s] failed to load the FM-index\n", __func__);
		return 1;
	}
	sr = rb3_srindex_build(&f, s_param, n_threads);
	if (sr == 0) {
		fprintf(stderr, "[E::%s] failed to build SR-index\n", __func__);
		rb3_fmi_free(&f);
		return 1;
	}
	if (fn == 0) {
		fn = RB3_MALLOC(char, strlen(argv[o.ind]) + 5);
		strcat(strcpy(fn, argv[o.ind]), ".sri");
		rb3_srindex_dump(sr, fn);
		free(fn);
	} else {
		rb3_srindex_dump(sr, fn);
	}
	if (rb3_verbose >= 3)
		fprintf(stderr, "[M::%s] SR-index built: %lld runs, s=%d, %lld subsampled, %lld sentinels\n",
		        __func__, (long long)sr->n_runs, sr->s, (long long)sr->n_sub, (long long)sr->m);
	rb3_fmi_free(&f);
	rb3_srindex_destroy(sr);
	return 0;
}
