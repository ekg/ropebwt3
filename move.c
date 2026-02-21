#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "rb3priv.h"
#include "move.h"
#include "lcp.h"
#include "rle.h"
#include "kalloc.h"

#define RB3_MVI_MAGIC "MVI\1"
#define RB3_MVI_HDR_SIZE 96
#define RB3_MVI_ROW_SIZE 48

static void move_push_run(rb3_move_t *m, int64_t *n_alloc, int64_t *cnt, int c, int64_t len, int64_t start)
{
	RB3_GROW(rb3_move_row_t, m->rows, m->n_runs, *n_alloc);
	rb3_move_row_t *row = &m->rows[m->n_runs++];
	memset(row, 0, sizeof(*row));
	row->c = c;
	row->len = len;
	row->p = start;
	row->pi = m->acc[c] + cnt[c];
	cnt[c] += len;
}

rb3_move_t *rb3_move_build(const rb3_fmi_t *f)
{
	rb3_move_t *m;
	int64_t i, pos, cnt[RB3_ASIZE], n_alloc = 0;
	int last_c = -1;
	int64_t run_start = 0, run_len = 0;

	m = RB3_CALLOC(rb3_move_t, 1);
	rb3_fmi_get_acc(f, m->acc);
	m->bwt_len = m->acc[RB3_ASIZE];
	m->n_runs = 0;
	m->rows = 0;
	m->d = 0;
	memset(cnt, 0, sizeof(cnt));
	pos = 0;

	/* Scan BWT and extract maximal runs. Adjacent same-character entries
	 * are merged (needed for FMR backend where runs may span blocks). */
	if (f->e) { /* FMD backend */
		rlditr_t itr;
		int c;
		int64_t l;
		rld_itr_init(f->e, &itr, 0);
		while ((l = rld_dec(f->e, &itr, &c, 0)) > 0) {
			if (c != last_c) {
				if (last_c >= 0)
					move_push_run(m, &n_alloc, cnt, last_c, run_len, run_start);
				last_c = c;
				run_start = pos;
				run_len = l;
			} else run_len += l;
			pos += l;
		}
	} else if (f->r) { /* FMR backend */
		mritr_t ri;
		const uint8_t *block;
		mr_itr_first(f->r, &ri, 0);
		while ((block = mr_itr_next_block(&ri)) != 0) {
			const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
			while (q < end) {
				int c = 0;
				int64_t l;
				rle_dec1(q, c, l);
				if (c != last_c) {
					if (last_c >= 0)
						move_push_run(m, &n_alloc, cnt, last_c, run_len, run_start);
					last_c = c;
					run_start = pos;
					run_len = l;
				} else run_len += l;
				pos += l;
			}
		}
	}
	if (last_c >= 0) /* flush the last run */
		move_push_run(m, &n_alloc, cnt, last_c, run_len, run_start);
	assert(pos == m->bwt_len);

	/* Compute destination run index xi: for each row i, find the row j
	 * whose run contains position pi (binary search on run start offsets). */
	for (i = 0; i < m->n_runs; ++i) {
		int64_t pi = m->rows[i].pi;
		int64_t lo = 0, hi = m->n_runs - 1;
		while (lo < hi) {
			int64_t mid = lo + (hi - lo + 1) / 2;
			if (m->rows[mid].p <= pi) lo = mid;
			else hi = mid - 1;
		}
		m->rows[i].xi = lo;
	}

	return m;
}

void rb3_move_destroy(rb3_move_t *m)
{
	if (m == 0) return;
	if (m->_mmap_base) {
		munmap(m->_mmap_base, m->_mmap_len);
		close(m->_fd);
	} else {
		free(m->rows);
	}
	free(m);
}

/*
 * Run splitting: split runs longer than a threshold so that fast-forward
 * during LF-mapping is bounded by < 2d steps. A run of length L is split
 * into ceil(L / max_len) subruns where max_len is derived from the splitting
 * depth d. After splitting, xi pointers are recomputed.
 */
void rb3_move_split(rb3_move_t *m, int32_t d)
{
	int64_t i, j, new_n, n_alloc;
	int64_t max_len;
	rb3_move_row_t *new_rows;

	if (d <= 0) return;
	m->d = d;

	/* Compute max run length: for depth d, threshold is n_runs^((d-1)/d).
	 * A simpler practical bound: split so each subrun has length <= max_len,
	 * where max_len = max(1, n_runs^((d-1)/d)). For small r, just use r. */
	{
		double base = (double)m->n_runs;
		double exp = (double)(d - 1) / (double)d;
		double v = 1.0;
		int k;
		for (k = 0; k < d - 1; k++) v *= base;
		/* v = base^(d-1), we want base^((d-1)/d) = v^(1/d) */
		/* Use pow approximation via repeated sqrt */
		if (d == 1) max_len = 1;
		else if (d == 2) {
			/* base^(1/2) = sqrt(base) */
			double s = base;
			/* Newton's method for integer sqrt */
			max_len = 1;
			while (max_len * max_len < (int64_t)s) max_len++;
		} else {
			/* For d=3: base^(2/3) */
			/* General: use exp/log approach */
			double logv = 0.0;
			double tmp = base;
			int iters = 0;
			/* Simple log2 via repeated halving */
			while (tmp > 2.0 && iters < 100) { tmp /= 2.0; logv += 1.0; iters++; }
			logv += (tmp - 1.0); /* rough log2 for remainder */
			logv *= exp;
			/* 2^logv */
			max_len = 1;
			while (logv >= 1.0) { max_len *= 2; logv -= 1.0; }
			if (logv > 0.0) max_len = (int64_t)(max_len * (1.0 + logv * 0.6931)); /* rough */
		}
		if (max_len < 1) max_len = 1;
	}

	/* First pass: count new rows */
	new_n = 0;
	for (i = 0; i < m->n_runs; i++) {
		int64_t L = m->rows[i].len;
		new_n += (L + max_len - 1) / max_len;
	}
	if (new_n == m->n_runs) return; /* nothing to split */

	/* Second pass: build new rows array */
	n_alloc = new_n;
	new_rows = RB3_MALLOC(rb3_move_row_t, n_alloc);
	j = 0;
	for (i = 0; i < m->n_runs; i++) {
		rb3_move_row_t *r = &m->rows[i];
		int64_t L = r->len;
		int64_t n_sub = (L + max_len - 1) / max_len;
		int64_t sub_len = L / n_sub; /* base length */
		int64_t extra = L - sub_len * n_sub; /* first 'extra' subruns get +1 */
		int64_t off_p = 0, off_pi = 0, k;
		for (k = 0; k < n_sub; k++) {
			int64_t slen = sub_len + (k < extra ? 1 : 0);
			memset(&new_rows[j], 0, sizeof(new_rows[j]));
			new_rows[j].c = r->c;
			new_rows[j].len = slen;
			new_rows[j].p = r->p + off_p;
			new_rows[j].pi = r->pi + off_pi;
			j++;
			off_p += slen;
			off_pi += slen;
		}
	}
	assert(j == new_n);

	free(m->rows);
	m->rows = new_rows;
	m->n_runs = new_n;

	/* Recompute xi via binary search */
	for (i = 0; i < m->n_runs; i++) {
		int64_t pi = m->rows[i].pi;
		int64_t lo = 0, hi = m->n_runs - 1;
		while (lo < hi) {
			int64_t mid = lo + (hi - lo + 1) / 2;
			if (m->rows[mid].p <= pi) lo = mid;
			else hi = mid - 1;
		}
		m->rows[i].xi = lo;
	}
}

/*
 * Precompute reposition distances: for each row i and each character c,
 * store the signed distance (in rows) to the nearest run of character c.
 * dist[c] = 0 means row i itself has character c.
 * dist[c] > 0 means go forward dist[c] rows.
 * dist[c] < 0 means go backward |dist[c]| rows.
 * If character c doesn't appear, dist[c] = 0 (will never be queried).
 */
void rb3_move_precompute_dist(rb3_move_t *m)
{
	int64_t i;
	int c;
	int64_t last_seen[RB3_ASIZE]; /* last row index where character c was seen */

	/* Forward pass: for each character, record distance to nearest previous occurrence */
	for (c = 0; c < RB3_ASIZE; c++)
		last_seen[c] = -1;
	for (i = 0; i < m->n_runs; i++) {
		int rc = m->rows[i].c;
		last_seen[rc] = i;
		for (c = 0; c < RB3_ASIZE; c++) {
			if (last_seen[c] >= 0)
				m->rows[i].dist[c] = (int16_t)(last_seen[c] - i); /* negative or zero */
			else
				m->rows[i].dist[c] = INT16_MAX; /* sentinel: not yet seen */
		}
	}

	/* Backward pass: check if a closer occurrence exists ahead */
	for (c = 0; c < RB3_ASIZE; c++)
		last_seen[c] = -1;
	for (i = m->n_runs - 1; i >= 0; i--) {
		int rc = m->rows[i].c;
		last_seen[rc] = i;
		for (c = 0; c < RB3_ASIZE; c++) {
			if (last_seen[c] < 0) continue; /* no occurrence ahead either */
			int64_t fwd_dist = last_seen[c] - i; /* positive or zero */
			int64_t cur_dist = m->rows[i].dist[c];
			if (cur_dist == INT16_MAX) {
				/* No backward occurrence, use forward */
				m->rows[i].dist[c] = (int16_t)fwd_dist;
			} else {
				/* cur_dist is <= 0 (backward). Compare |cur_dist| vs fwd_dist */
				if (fwd_dist < -cur_dist)
					m->rows[i].dist[c] = (int16_t)fwd_dist;
				/* else keep backward (ties go to backward, which is already stored) */
			}
		}
	}

	/* Any remaining INT16_MAX means character doesn't exist at all; set to 0 */
	for (i = 0; i < m->n_runs; i++)
		for (c = 0; c < RB3_ASIZE; c++)
			if (m->rows[i].dist[c] == INT16_MAX)
				m->rows[i].dist[c] = 0;
}

/*
 * LF-mapping via move table.
 * Given BWT position pos in run *run_idx, compute LF(pos).
 * The formula is: LF(pos) = M[i].pi + (pos - M[i].p)
 * Then follow xi to the destination run and fast-forward to find
 * the exact run containing the result.
 */
int64_t rb3_move_lf(const rb3_move_t *m, int64_t pos, int64_t *run_idx)
{
	int64_t i = *run_idx;
	int64_t lf_pos;
	int64_t dest;

	/* Compute LF(pos) */
	lf_pos = m->rows[i].pi + (pos - m->rows[i].p);

	/* Jump to destination row */
	dest = m->rows[i].xi;

	/* Fast-forward: scan forward to find the run containing lf_pos */
	while (dest + 1 < m->n_runs && m->rows[dest + 1].p <= lf_pos)
		dest++;
	/* Also check backward (shouldn't happen with correct xi, but be safe) */
	while (dest > 0 && m->rows[dest].p > lf_pos)
		dest--;

	*run_idx = dest;
	return lf_pos;
}

/*
 * Reposition: when the current run's character doesn't match the query
 * character c, jump to the nearest run of character c using precomputed
 * distances.
 * Returns the new run index. The caller should use rows[new_idx].p as
 * the new BWT position (head of the target run).
 */
int64_t rb3_move_reposition(const rb3_move_t *m, int64_t run_idx, int8_t c)
{
	return run_idx + m->rows[run_idx].dist[c];
}

/*
 * Combined backward search step: reposition to character c, then LF-map.
 * This is the core primitive for backward search on the move structure.
 *
 * Given current BWT position pos in run *run_idx, and query character c:
 * 1. If the current run's character != c, reposition to nearest run of c
 * 2. Apply LF-mapping from the (possibly repositioned) position
 *
 * Returns the new BWT position and updates *run_idx.
 */
int64_t rb3_move_step(const rb3_move_t *m, int64_t pos, int64_t *run_idx, int8_t c)
{
	int64_t i = *run_idx;

	/* Reposition if character doesn't match */
	if (m->rows[i].c != c) {
		i = rb3_move_reposition(m, i, c);
		pos = m->rows[i].p; /* move to head of the target run */
	}

	*run_idx = i;
	return rb3_move_lf(m, pos, run_idx);
}

/*
 * Binary search for the run containing BWT position pos.
 * Returns the index i such that rows[i].p <= pos < rows[i].p + rows[i].len,
 * or equivalently the largest i with rows[i].p <= pos.
 */
static int64_t move_find_run(const rb3_move_t *m, int64_t pos)
{
	int64_t lo = 0, hi = m->n_runs - 1;
	while (lo < hi) {
		int64_t mid = lo + (hi - lo + 1) / 2;
		if (m->rows[mid].p <= pos) lo = mid;
		else hi = mid - 1;
	}
	return lo;
}

/*
 * Count occurrences of pattern[0..len-1] using interval-tracking backward search.
 * Pattern is encoded as 0-5 integers (same as BWT alphabet).
 *
 * Uses a precomputed cumulative rank table for O(1) rank queries per step.
 * Total time: O(r) setup + O(|P| * log r) search.
 */
int64_t rb3_move_count(const rb3_move_t *m, int len, const uint8_t *pattern)
{
	int64_t lo, hi, *cumrank;
	int i;
	int8_t c;

	if (len <= 0) return m->bwt_len;
	if (m->n_runs == 0) return 0;

	/* Build cumulative rank table: cumrank[i * 6 + c] = rank(c, rows[i].p)
	 * i.e., the number of character c in BWT[0..rows[i].p). */
	cumrank = RB3_CALLOC(int64_t, ((int64_t)m->n_runs + 1) * RB3_ASIZE);
	{
		int64_t j;
		for (j = 0; j < m->n_runs; j++) {
			int c2;
			for (c2 = 0; c2 < RB3_ASIZE; c2++)
				cumrank[(j + 1) * RB3_ASIZE + c2] = cumrank[j * RB3_ASIZE + c2];
			cumrank[(j + 1) * RB3_ASIZE + m->rows[j].c] += m->rows[j].len;
		}
	}

	/* Initialize with last character of pattern */
	c = pattern[len - 1];
	if (c < 0 || c >= RB3_ASIZE || m->acc[c] >= m->acc[c + 1]) {
		free(cumrank);
		return 0;
	}
	lo = m->acc[c];
	hi = m->acc[c + 1];

	/* Extend backward through the pattern */
	for (i = len - 2; i >= 0; --i) {
		int64_t lo_run, hi_run, rank_lo, rank_hi;
		c = pattern[i];
		if (c < 0 || c >= RB3_ASIZE || m->acc[c] >= m->acc[c + 1]) {
			free(cumrank);
			return 0;
		}

		/* Find runs containing lo and hi */
		lo_run = move_find_run(m, lo);
		hi_run = (hi < m->bwt_len) ? move_find_run(m, hi) : m->n_runs - 1;

		/* rank(c, pos) = cumrank[run][c] + (run's char == c ? pos - run.p : 0) */
		rank_lo = cumrank[lo_run * RB3_ASIZE + c];
		if (m->rows[lo_run].c == c)
			rank_lo += lo - m->rows[lo_run].p;

		if (hi >= m->bwt_len) {
			rank_hi = cumrank[m->n_runs * RB3_ASIZE + c];
		} else {
			rank_hi = cumrank[hi_run * RB3_ASIZE + c];
			if (m->rows[hi_run].c == c)
				rank_hi += hi - m->rows[hi_run].p;
		}

		lo = m->acc[c] + rank_lo;
		hi = m->acc[c] + rank_hi;
		if (lo >= hi) {
			free(cumrank);
			return 0;
		}
	}

	free(cumrank);
	return hi - lo;
}

/*
 * Checksum: XOR all 64-bit words in the row data.
 */
static uint64_t mvi_checksum(const rb3_move_row_t *rows, int64_t n_runs)
{
	const uint64_t *p = (const uint64_t *)rows;
	int64_t n_words = n_runs * (int64_t)(sizeof(rb3_move_row_t) / sizeof(uint64_t));
	uint64_t cs = 0;
	int64_t i;
	for (i = 0; i < n_words; i++)
		cs ^= p[i];
	return cs;
}

/*
 * Save move table to a .mvi binary file.
 *
 * Header layout (96 bytes):
 *   [0:4]   magic "MVI\1"
 *   [4:8]   uint32_t flags (reserved)
 *   [8:16]  int64_t n_runs
 *   [16:24] int64_t bwt_len
 *   [24:80] int64_t acc[7]
 *   [80:84] int32_t d
 *   [84:88] uint32_t row_size
 *   [88:96] uint64_t checksum
 *
 * Body: n_runs × rb3_move_row_t (48 bytes each)
 */
int rb3_move_save(const rb3_move_t *m, const char *fn)
{
	FILE *fp;
	uint32_t flags = 0, row_size = RB3_MVI_ROW_SIZE;
	uint64_t checksum;

	assert(sizeof(rb3_move_row_t) == RB3_MVI_ROW_SIZE);

	fp = fopen(fn, "wb");
	if (fp == 0) return -1;

	fwrite(RB3_MVI_MAGIC, 1, 4, fp);
	fwrite(&flags, 4, 1, fp);
	fwrite(&m->n_runs, 8, 1, fp);
	fwrite(&m->bwt_len, 8, 1, fp);
	fwrite(m->acc, 8, 7, fp);
	fwrite(&m->d, 4, 1, fp);
	fwrite(&row_size, 4, 1, fp);
	checksum = mvi_checksum(m->rows, m->n_runs);
	fwrite(&checksum, 8, 1, fp);

	fwrite(m->rows, sizeof(rb3_move_row_t), m->n_runs, fp);
	fclose(fp);
	return 0;
}

/*
 * Load move table from a .mvi file via memory mapping.
 * Returns NULL on error. The loaded rows are read-only (mmap'd).
 */
rb3_move_t *rb3_move_load(const char *fn)
{
	rb3_move_t *m;
	int fd;
	struct stat st;
	uint8_t *base;
	uint32_t flags, row_size;
	uint64_t checksum, computed_cs;
	size_t expected_size;

	fd = open(fn, O_RDONLY);
	if (fd < 0) return 0;
	if (fstat(fd, &st) != 0 || st.st_size < RB3_MVI_HDR_SIZE) { close(fd); return 0; }

	base = (uint8_t *)mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (base == MAP_FAILED) { close(fd); return 0; }

	/* Validate magic */
	if (memcmp(base, RB3_MVI_MAGIC, 4) != 0) goto fail;

	m = RB3_CALLOC(rb3_move_t, 1);

	/* Parse header */
	memcpy(&flags, base + 4, 4);
	memcpy(&m->n_runs, base + 8, 8);
	memcpy(&m->bwt_len, base + 16, 8);
	memcpy(m->acc, base + 24, 56);
	memcpy(&m->d, base + 80, 4);
	memcpy(&row_size, base + 84, 4);
	memcpy(&checksum, base + 88, 8);

	if (row_size != RB3_MVI_ROW_SIZE) goto fail2;

	expected_size = RB3_MVI_HDR_SIZE + (size_t)m->n_runs * row_size;
	if ((size_t)st.st_size != expected_size) goto fail2;

	m->rows = (rb3_move_row_t *)(base + RB3_MVI_HDR_SIZE);

	/* Verify checksum */
	computed_cs = mvi_checksum(m->rows, m->n_runs);
	if (computed_cs != checksum) goto fail2;

	m->_fd = fd;
	m->_mmap_base = base;
	m->_mmap_len = st.st_size;
	return m;

fail2:
	free(m);
fail:
	munmap(base, st.st_size);
	close(fd);
	return 0;
}

/***********************************************
 * Move + LCP matching statistics (MONI-style) *
 ***********************************************/

/* rb3_bmove_rank1a is defined later in this file (public, declared in move.h) */

/* Find the LCP run index containing BWT position pos */
static inline int64_t move_lcp_find_run(const rb3_lcp_t *lcp, int64_t pos)
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
 * Precompute per-move-row thresholds from LCP thresholds.
 * Each move row corresponds to a sub-run (or whole run) of the BWT.
 * We map each move row to the LCP run containing its BWT position via
 * a linear scan (both arrays are sorted by BWT offset).
 * Returns a malloc'd array of m->n_runs int64_t values.
 */
int64_t *rb3_move_lcp_thresholds(const rb3_move_t *m, const rb3_lcp_t *lcp)
{
	int64_t *th, i, j;
	if (m == 0 || lcp == 0 || lcp->thresholds == 0) return 0;
	th = RB3_MALLOC(int64_t, m->n_runs);
	j = 0;
	for (i = 0; i < m->n_runs; ++i) {
		while (j + 1 < lcp->n_runs && lcp->run_starts[j + 1] <= m->rows[i].p)
			++j;
		th[i] = lcp->thresholds[j];
	}
	return th;
}

/*
 * Map each move row to its LCP run index.
 * Returns a malloc'd array of m->n_runs int64_t values.
 */
int64_t *rb3_move_lcp_run_map(const rb3_move_t *m, const rb3_lcp_t *lcp)
{
	int64_t *rm, i, j;
	if (m == 0 || lcp == 0) return 0;
	rm = RB3_MALLOC(int64_t, m->n_runs);
	j = 0;
	for (i = 0; i < m->n_runs; ++i) {
		while (j + 1 < lcp->n_runs && lcp->run_starts[j + 1] <= m->rows[i].p)
			++j;
		rm[i] = j;
	}
	return rm;
}

/*
 * One backward step of matching statistics using move + LCP.
 *
 * Given the current BWT position pos in run *run_idx with current match
 * length *match_len, extend backward with character c:
 *
 *   - If the current run's character equals c (match):
 *     Apply LF-mapping. match_len += 1.
 *
 *   - If the current run's character differs (mismatch):
 *     Reposition to the nearest run of character c. Compute the
 *     threshold as the range minimum of lcp_samples over all run
 *     boundaries crossed during reposition. Truncate match_len to
 *     min(match_len, threshold), then apply LF-mapping and increment.
 *
 * Returns the new BWT position, updates *run_idx and *match_len.
 * Returns -1 if character c does not exist in the BWT.
 */
int64_t rb3_move_ms_step(const rb3_move_t *m, const int64_t *run_map,
                         const rb3_lcp_t *lcp,
                         int64_t pos, int64_t *run_idx, int64_t *match_len, int8_t c)
{
	int64_t i = *run_idx;
	if (c < 1 || c >= RB3_ASIZE || m->acc[c] == m->acc[c + 1])
		return -1;
	if (m->rows[i].c == c) {
		/* Match: standard LF-mapping */
		pos = rb3_move_lf(m, pos, run_idx);
		++(*match_len);
	} else {
		/* Mismatch: reposition, compute threshold, truncate, LF-map.
		 *
		 * Direction-dependent threshold using MONI tau:
		 *   - If zone matches reposition direction (left zone going left,
		 *     right zone going right), the boundary LCP is exact.
		 *   - If zone/direction mismatch, fall back to within_min (the
		 *     minimum of all within-run LCPs, always a valid lower bound).
		 * Then combine with range-min of boundary LCPs crossed. */
		int64_t old_lcp_run, new_lcp_run, lo, hi, j, threshold;
		old_lcp_run = run_map[i];
		i = rb3_move_reposition(m, i, c);
		new_lcp_run = run_map[i];
		/* Direction-dependent threshold from tau */
		if (new_lcp_run < old_lcp_run) {
			/* Going LEFT */
			if (pos < lcp->tau[old_lcp_run])
				threshold = lcp->lcp_samples[old_lcp_run]; /* left zone: exact */
			else
				threshold = lcp->within_min[old_lcp_run]; /* right zone going left: fallback */
		} else {
			/* Going RIGHT (or same run after split) */
			if (pos >= lcp->tau[old_lcp_run])
				threshold = (old_lcp_run + 1 < lcp->n_runs) ? lcp->lcp_samples[old_lcp_run + 1] : 0; /* right zone: exact */
			else
				threshold = lcp->within_min[old_lcp_run]; /* left zone going right: fallback */
		}
		/* Range minimum of lcp_samples over boundaries between old and new runs */
		lo = old_lcp_run < new_lcp_run ? old_lcp_run : new_lcp_run;
		hi = old_lcp_run > new_lcp_run ? old_lcp_run : new_lcp_run;
		for (j = lo + 1; j <= hi; ++j)
			if (lcp->lcp_samples[j] < threshold)
				threshold = lcp->lcp_samples[j];
		if (*match_len > threshold)
			*match_len = threshold;
		pos = m->rows[i].p;
		*run_idx = i;
		pos = rb3_move_lf(m, pos, run_idx);
		++(*match_len);
	}
	return pos;
}

/*
 * Compute matching statistics for pattern[0..len-1] using move + LCP.
 *
 * MS[i] = length of the longest prefix of pattern[i..len-1] that occurs
 * as a substring of the reference text encoded in the BWT.
 *
 * Uses the same interval-based algorithm as the BWT version (rb3_ms_compute)
 * but replaces direct BWT rank queries with cumulative-rank lookups on the
 * b-move structure. This gives exact results.
 *
 * Returns 0 on success, -1 on error.
 */
int rb3_move_ms_compute(const rb3_move_t *m, const rb3_lcp_t *lcp,
                        int64_t len, const uint8_t *pattern, int64_t *ms)
{
	rb3_bmove_t *bm;
	int64_t i, k, l, d;
	int64_t ok[RB3_ASIZE], ol[RB3_ASIZE];

	if (len <= 0) return 0;
	if (m == 0 || m->n_runs == 0) return -1;
	if (lcp == 0 || lcp->thresholds == 0) return -1;

	bm = rb3_bmove_init(m);
	if (bm == 0) return -1;

	k = 0;
	l = m->bwt_len;
	d = 0;

	for (i = len - 1; i >= 0; --i) {
		int c = pattern[i];
		int64_t nk, nl;

		rb3_bmove_rank1a(bm, k, ok);
		rb3_bmove_rank1a(bm, l, ol);
		nk = m->acc[c] + ok[c];
		nl = m->acc[c] + ol[c];

		if (nk < nl) {
			k = nk; l = nl; ++d;
		} else {
			while (d > 0) {
				int64_t lcp_k, lcp_l, th;
				int64_t run_idx, lo_run, hi_run;
				int fc;

				lcp_k = (k > 0) ? rb3_lcp_at_position(lcp, k) : 0;
				lcp_l = (l > 0 && l < m->bwt_len) ? rb3_lcp_at_position(lcp, l) : 0;
				th = lcp_k > lcp_l ? lcp_k : lcp_l;

				if (th < d)
					d = th;
				else
					--d;

				for (fc = 0; fc < RB3_ASIZE; ++fc)
					if (k < m->acc[fc + 1]) break;

				run_idx = move_lcp_find_run(lcp, k);
				lo_run = run_idx;
				hi_run = move_lcp_find_run(lcp, l > 0 ? l - 1 : 0);
				while (lo_run > 0 && lcp->lcp_samples[lo_run] >= d)
					--lo_run;
				while (hi_run + 1 < lcp->n_runs && lcp->lcp_samples[hi_run + 1] >= d)
					++hi_run;
				k = lcp->run_starts[lo_run];
				if (hi_run + 1 < lcp->n_runs)
					l = lcp->run_starts[hi_run + 1];
				else
					l = m->bwt_len;

				if (d > 0) {
					if (k < m->acc[fc]) k = m->acc[fc];
					if (l > m->acc[fc + 1]) l = m->acc[fc + 1];
				}

				if (d == 0) break;

				rb3_bmove_rank1a(bm, k, ok);
				rb3_bmove_rank1a(bm, l, ol);
				nk = m->acc[c] + ok[c];
				nl = m->acc[c] + ol[c];
				if (nk < nl) { k = nk; l = nl; ++d; break; }
			}
			if (d == 0) {
				k = m->acc[c]; l = m->acc[c + 1];
				if (k < l) d = 1;
			}
		}
		ms[i] = d;
	}

	rb3_bmove_destroy(bm);
	return 0;
}

/***************************************************
 * b-move: bidirectional extension via move structure
 ***************************************************/

rb3_bmove_t *rb3_bmove_init(const rb3_move_t *mv)
{
	rb3_bmove_t *bm;
	int64_t i;
	int c;

	if (mv == 0 || mv->n_runs == 0) return 0;

	bm = RB3_CALLOC(rb3_bmove_t, 1);
	bm->mv = mv;
	bm->cumrank = RB3_CALLOC(int64_t, ((int64_t)mv->n_runs + 1) * RB3_ASIZE);

	for (i = 0; i < mv->n_runs; i++) {
		for (c = 0; c < RB3_ASIZE; c++)
			bm->cumrank[(i + 1) * RB3_ASIZE + c] = bm->cumrank[i * RB3_ASIZE + c];
		bm->cumrank[(i + 1) * RB3_ASIZE + mv->rows[i].c] += mv->rows[i].len;
	}

	return bm;
}

void rb3_bmove_destroy(rb3_bmove_t *bm)
{
	if (bm == 0) return;
	free(bm->cumrank);
	free(bm);
}

/*
 * Compute rank(c, pos) for all characters c using the cumulative rank table.
 * ok[c] = number of character c in BWT[0..pos).
 */
void rb3_bmove_rank1a(const rb3_bmove_t *bm, int64_t pos, int64_t ok[RB3_ASIZE])
{
	const rb3_move_t *m = bm->mv;
	int c;
	int64_t run;

	if (pos <= 0) {
		memset(ok, 0, RB3_ASIZE * sizeof(int64_t));
		return;
	}
	if (pos >= m->bwt_len) {
		for (c = 0; c < RB3_ASIZE; c++)
			ok[c] = bm->cumrank[m->n_runs * RB3_ASIZE + c];
		return;
	}

	run = move_find_run(m, pos);
	for (c = 0; c < RB3_ASIZE; c++)
		ok[c] = bm->cumrank[run * RB3_ASIZE + c];
	ok[m->rows[run].c] += pos - m->rows[run].p;
}

/*
 * Dual rank query: compute rank arrays at positions k and l.
 */
void rb3_bmove_rank2a(const rb3_bmove_t *bm, int64_t k, int64_t l, int64_t ok[RB3_ASIZE], int64_t ol[RB3_ASIZE])
{
	rb3_bmove_rank1a(bm, k, ok);
	rb3_bmove_rank1a(bm, l, ol);
}

/*
 * FMD-style bidirectional extension using b-move.
 *
 * Replaces rb3_fmd_extend(): instead of rank queries on the FMD/FMR BWT,
 * uses the cumulative rank table for O(log r) rank at arbitrary positions.
 * The FMD symmetry (both DNA strands in one BWT) is preserved.
 */
void rb3_bmove_extend(const rb3_bmove_t *bm, const rb3_sai_t *ik, rb3_sai_t ok[RB3_ASIZE], int is_back)
{
	int64_t tk[RB3_ASIZE], tl[RB3_ASIZE];
	int c;

	is_back = !!is_back;
	rb3_bmove_rank1a(bm, ik->x[!is_back], tk);
	rb3_bmove_rank1a(bm, ik->x[!is_back] + ik->size, tl);

	for (c = 0; c < RB3_ASIZE; c++) {
		ok[c].x[!is_back] = bm->mv->acc[c] + tk[c];
		ok[c].size = (tl[c] -= tk[c]);
	}
	ok[0].x[is_back] = ik->x[is_back];
	ok[4].x[is_back] = ok[0].x[is_back] + tl[0];
	ok[3].x[is_back] = ok[4].x[is_back] + tl[4];
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[5].x[is_back] = ok[1].x[is_back] + tl[1];
}

/*
 * Helper: initialize a bidirectional interval for a single character.
 * Same as rb3_fmd_set_intv() but uses the move structure's acc[].
 */
static inline void bmove_set_intv(const rb3_bmove_t *bm, int c, rb3_sai_t *ik)
{
	const int64_t *acc = bm->mv->acc;
	ik->x[0] = acc[c];
	ik->size = acc[c + 1] - acc[c];
	ik->x[1] = acc[rb3_comp(c)];
	ik->info = 0;
}

static void bmove_sai_reverse(rb3_sai_t *a, int64_t l)
{
	int64_t i;
	rb3_sai_t t;
	for (i = 0; i < l>>1; ++i)
		t = a[i], a[i] = a[l - 1 - i], a[l - 1 - i] = t;
}

/*
 * SMEM finding for one seed position using b-move (original algorithm).
 * Same logic as rb3_fmd_smem1() but uses rb3_bmove_extend().
 */
static int64_t rb3_bmove_smem1(void *km, const rb3_bmove_t *bm, int64_t min_occ, int64_t min_len, int64_t len, const uint8_t *q, int64_t x, rb3_sai_v *mem, rb3_sai_v *curr, rb3_sai_v *prev)
{
	int64_t i, j, ret;
	rb3_sai_t ik, ok[6];
	rb3_sai_v *swap;
	size_t oldn = mem->n;

	assert(len <= INT32_MAX);
	bmove_set_intv(bm, q[x], &ik);
	ik.info = x + 1;
	if (ik.size == 0) return x + 1;
	for (i = x + 1, curr->n = 0; i < len; ++i) { /* forward extension */
		int c = rb3_comp(q[i]);
		rb3_bmove_extend(bm, &ik, ok, 0);
		if (ok[c].size != ik.size) {
			Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
			curr->a[curr->n++] = ik;
			if (ok[c].size < min_occ) break;
		}
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) {
		Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
		curr->a[curr->n++] = ik;
	}
	bmove_sai_reverse(curr->a, curr->n);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { /* backward extension */
		int c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < (int64_t)prev->n; ++j) {
			rb3_sai_t *p = &prev->a[j];
			rb3_bmove_extend(bm, p, ok, 1);
			if (c == 0 || ok[c].size < min_occ) {
				if (curr->n == 0 && (int32_t)p->info - i - 1 >= min_len && (mem->n == oldn || i + 1 < mem->a[mem->n-1].info>>32)) {
					rb3_sai_t *q2;
					Kgrow(km, rb3_sai_t, mem->a, mem->n, mem->m);
					q2 = &mem->a[mem->n++];
					*q2 = *p; q2->info |= (int64_t)(i + 1)<<32;
				}
			} else if (curr->n == 0 || ok[c].size != curr->a[curr->n-1].size) {
				ok[c].info = p->info;
				Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
				curr->a[curr->n++] = ok[c];
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	bmove_sai_reverse(&mem->a[oldn], mem->n - oldn);
	return ret;
}

int64_t rb3_bmove_smem(void *km, const rb3_bmove_t *bm, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len)
{
	int64_t x = 0;
	rb3_sai_v curr = {0,0,0}, prev = {0,0,0};
	mem->n = 0;
	do {
		x = rb3_bmove_smem1(km, bm, min_occ, min_len, len, q, x, mem, &curr, &prev);
	} while (x < len);
	kfree(km, curr.a);
	kfree(km, prev.a);
	return mem->n;
}

/*
 * SMEM finding for one seed position using b-move (Travis Gagie algorithm).
 * Same logic as rb3_fmd_smem1_TG() but uses rb3_bmove_extend().
 */
static int64_t rb3_bmove_smem1_TG(void *km, const rb3_bmove_t *bm, int64_t min_occ, int64_t min_len, int64_t len, const uint8_t *q, int64_t x, rb3_sai_v *mem, int32_t check_long)
{
	int64_t i, j;
	rb3_sai_t ik, ok[6];

	assert(len <= INT32_MAX);
	if (len - x < min_len) return len;
	bmove_set_intv(bm, q[x + min_len - 1], &ik);
	for (i = x + min_len - 2; i >= x; --i) { /* backward extension */
		int c = q[i];
		rb3_bmove_extend(bm, &ik, ok, 1);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	if (i >= x) return i + 1; /* no MEM found */
	if (check_long) return -1;
	for (j = x + min_len; j < len; ++j) { /* forward extension */
		int c = rb3_comp(q[j]);
		rb3_bmove_extend(bm, &ik, ok, 0);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	{
		rb3_sai_t *p;
		Kgrow(km, rb3_sai_t, mem->a, mem->n, mem->m);
		p = &mem->a[mem->n++];
		*p = ik;
		p->info = (uint64_t)x<<32 | j;
	}
	if (j == len) return len;
	bmove_set_intv(bm, q[j], &ik);
	for (i = j - 1; i > x; --i) { /* backward extension again */
		int c = q[i];
		rb3_bmove_extend(bm, &ik, ok, 1);
		if (ok[c].size < min_occ) break;
		ik = ok[c];
	}
	return i + 1;
}

int64_t rb3_bmove_smem_TG(void *km, const rb3_bmove_t *bm, int64_t len, const uint8_t *q, rb3_sai_v *mem, int64_t min_occ, int64_t min_len)
{
	int64_t x = 0;
	mem->n = 0;
	do {
		x = rb3_bmove_smem1_TG(km, bm, min_occ, min_len, len, q, x, mem, 0);
	} while (x < len);
	return mem->n;
}
