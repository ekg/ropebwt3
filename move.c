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

/*
 * .mvi file format
 *
 * V1 (legacy): magic "MVI\1", 96-byte header, n_runs × 48-byte rows
 * V2 (compact): magic "MVI\2", 96-byte header, then:
 *     uint32_t xi[n_runs]
 *     uint16_t len[n_runs]   (requires run splitting, all lens < 65536)
 *     int8_t   c[n_runs]
 *     int16_t  dist[n_runs * 6]
 *   p[] and pi[] are derived at load time from prefix sums.
 *   Total: 7 + 12 = 19 bytes/row on disk (2.5x smaller than v1).
 */
#define RB3_MVI_MAGIC_V1 "MVI\1"
#define RB3_MVI_MAGIC_V2 "MVI\2"
#define RB3_MVI_HDR_SIZE 96
#define RB3_MVI_ROW_SIZE_V1 48

/* Allocate separate arrays for n_runs rows */
static void move_alloc_arrays(rb3_move_t *m, int64_t n)
{
	m->p    = RB3_CALLOC(int64_t, n);
	m->pi   = RB3_CALLOC(int64_t, n);
	m->xi   = RB3_CALLOC(uint32_t, n);
	m->len  = RB3_CALLOC(int64_t, n);
	m->c    = RB3_CALLOC(int8_t, n);
	m->dist = RB3_CALLOC(int16_t, n * RB3_ASIZE);
}

static void move_free_arrays(rb3_move_t *m)
{
	free(m->p);
	free(m->pi);
	free(m->xi);
	free(m->len);
	free(m->c);
	free(m->dist);
}

/* Derive p[] from prefix sum of len[] */
static void move_derive_p(rb3_move_t *m)
{
	int64_t i, pos = 0;
	for (i = 0; i < m->n_runs; i++) {
		m->p[i] = pos;
		pos += m->len[i];
	}
}

/* Derive pi[] from acc[] and cumulative character counts */
static void move_derive_pi(rb3_move_t *m)
{
	int64_t i, cnt[RB3_ASIZE];
	memset(cnt, 0, sizeof(cnt));
	for (i = 0; i < m->n_runs; i++) {
		int cc = m->c[i];
		m->pi[i] = m->acc[cc] + cnt[cc];
		cnt[cc] += m->len[i];
	}
}

/* Compute xi[] via binary search on p[] */
static void move_compute_xi(rb3_move_t *m)
{
	int64_t i;
	for (i = 0; i < m->n_runs; i++) {
		int64_t target = m->pi[i];
		int64_t lo = 0, hi = m->n_runs - 1;
		while (lo < hi) {
			int64_t mid = lo + (hi - lo + 1) / 2;
			if (m->p[mid] <= target) lo = mid;
			else hi = mid - 1;
		}
		m->xi[i] = (uint32_t)lo;
	}
}

static void move_push_run(rb3_move_t *m, int64_t *n_alloc, int64_t *cnt, int cc, int64_t run_len, int64_t start)
{
	int64_t i = m->n_runs;
	if (i >= *n_alloc) {
		int64_t new_alloc = i + 1;
		new_alloc += (new_alloc >> 1) + 16;
		m->p   = RB3_REALLOC(int64_t, m->p, new_alloc);
		m->pi  = RB3_REALLOC(int64_t, m->pi, new_alloc);
		m->xi  = RB3_REALLOC(uint32_t, m->xi, new_alloc);
		m->len = RB3_REALLOC(int64_t, m->len, new_alloc);
		m->c   = RB3_REALLOC(int8_t, m->c, new_alloc);
		/* dist is allocated separately after build */
		*n_alloc = new_alloc;
	}
	m->p[i] = start;
	m->c[i] = cc;
	m->len[i] = run_len;
	m->pi[i] = m->acc[cc] + cnt[cc];
	m->xi[i] = 0; /* computed later */
	cnt[cc] += run_len;
	m->n_runs++;
}

rb3_move_t *rb3_move_build(const rb3_fmi_t *f)
{
	rb3_move_t *m;
	int64_t pos, cnt[RB3_ASIZE], n_alloc = 0;
	int last_c = -1;
	int64_t run_start = 0, run_len = 0;

	m = RB3_CALLOC(rb3_move_t, 1);
	rb3_fmi_get_acc(f, m->acc);
	m->bwt_len = m->acc[RB3_ASIZE];
	m->n_runs = 0;
	m->d = 0;
	memset(cnt, 0, sizeof(cnt));
	pos = 0;

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
	if (last_c >= 0)
		move_push_run(m, &n_alloc, cnt, last_c, run_len, run_start);
	assert(pos == m->bwt_len);

	/* Allocate dist array (initially zero) */
	m->dist = RB3_CALLOC(int16_t, m->n_runs * RB3_ASIZE);

	/* Compute destination run index xi via binary search on p[] */
	move_compute_xi(m);

	return m;
}

void rb3_move_destroy(rb3_move_t *m)
{
	if (m == 0) return;
	if (m->_mmap_base) {
		munmap(m->_mmap_base, m->_mmap_len);
		close(m->_fd);
	}
	move_free_arrays(m);
	free(m);
}

/*
 * Run splitting: split runs longer than a threshold so that fast-forward
 * during LF-mapping is bounded by < 2d steps.
 */
void rb3_move_split(rb3_move_t *m, int32_t d)
{
	int64_t i, j, new_n, max_len;
	int64_t *new_p, *new_pi, *new_len;
	uint32_t *new_xi;
	int8_t *new_c;

	if (d <= 0) return;
	m->d = d;

	/* Compute max run length */
	{
		double base = (double)m->n_runs;
		if (d == 1) max_len = 1;
		else if (d == 2) {
			max_len = 1;
			while (max_len * max_len < (int64_t)base) max_len++;
		} else {
			double logv = 0.0, tmp = base, exp = (double)(d - 1) / (double)d;
			int iters = 0;
			while (tmp > 2.0 && iters < 100) { tmp /= 2.0; logv += 1.0; iters++; }
			logv += (tmp - 1.0);
			logv *= exp;
			max_len = 1;
			while (logv >= 1.0) { max_len *= 2; logv -= 1.0; }
			if (logv > 0.0) max_len = (int64_t)(max_len * (1.0 + logv * 0.6931));
		}
		if (max_len < 1) max_len = 1;
	}

	/* First pass: count new rows */
	new_n = 0;
	for (i = 0; i < m->n_runs; i++)
		new_n += (m->len[i] + max_len - 1) / max_len;
	if (new_n == m->n_runs) return;

	/* Second pass: build new arrays */
	new_p   = RB3_MALLOC(int64_t, new_n);
	new_pi  = RB3_MALLOC(int64_t, new_n);
	new_xi  = RB3_MALLOC(uint32_t, new_n);
	new_len = RB3_MALLOC(int64_t, new_n);
	new_c   = RB3_MALLOC(int8_t, new_n);
	j = 0;
	for (i = 0; i < m->n_runs; i++) {
		int64_t L = m->len[i];
		int64_t n_sub = (L + max_len - 1) / max_len;
		int64_t sub_len = L / n_sub;
		int64_t extra = L - sub_len * n_sub;
		int64_t off_p = 0, off_pi = 0, k;
		for (k = 0; k < n_sub; k++) {
			int64_t slen = sub_len + (k < extra ? 1 : 0);
			new_c[j]   = m->c[i];
			new_len[j] = slen;
			new_p[j]   = m->p[i] + off_p;
			new_pi[j]  = m->pi[i] + off_pi;
			new_xi[j]  = 0; /* recomputed below */
			j++;
			off_p += slen;
			off_pi += slen;
		}
	}
	assert(j == new_n);

	free(m->p); free(m->pi); free(m->xi); free(m->len); free(m->c);
	m->p = new_p; m->pi = new_pi; m->xi = new_xi; m->len = new_len; m->c = new_c;
	m->n_runs = new_n;

	/* Reallocate dist (will be recomputed by caller) */
	free(m->dist);
	m->dist = RB3_CALLOC(int16_t, new_n * RB3_ASIZE);

	/* Recompute xi via binary search */
	move_compute_xi(m);
}

void rb3_move_precompute_dist(rb3_move_t *m)
{
	int64_t i;
	int cc;
	int64_t last_seen[RB3_ASIZE];

	/* Forward pass */
	for (cc = 0; cc < RB3_ASIZE; cc++)
		last_seen[cc] = -1;
	for (i = 0; i < m->n_runs; i++) {
		last_seen[(int)m->c[i]] = i;
		for (cc = 0; cc < RB3_ASIZE; cc++) {
			if (last_seen[cc] >= 0)
				m->dist[i * RB3_ASIZE + cc] = (int16_t)(last_seen[cc] - i);
			else
				m->dist[i * RB3_ASIZE + cc] = INT16_MAX;
		}
	}

	/* Backward pass */
	for (cc = 0; cc < RB3_ASIZE; cc++)
		last_seen[cc] = -1;
	for (i = m->n_runs - 1; i >= 0; i--) {
		last_seen[(int)m->c[i]] = i;
		for (cc = 0; cc < RB3_ASIZE; cc++) {
			if (last_seen[cc] < 0) continue;
			int64_t fwd_dist = last_seen[cc] - i;
			int64_t cur_dist = m->dist[i * RB3_ASIZE + cc];
			if (cur_dist == INT16_MAX) {
				m->dist[i * RB3_ASIZE + cc] = (int16_t)fwd_dist;
			} else {
				if (fwd_dist < -cur_dist)
					m->dist[i * RB3_ASIZE + cc] = (int16_t)fwd_dist;
			}
		}
	}

	for (i = 0; i < m->n_runs; i++)
		for (cc = 0; cc < RB3_ASIZE; cc++)
			if (m->dist[i * RB3_ASIZE + cc] == INT16_MAX)
				m->dist[i * RB3_ASIZE + cc] = 0;
}

int64_t rb3_move_lf(const rb3_move_t *m, int64_t pos, int64_t *run_idx)
{
	int64_t i = *run_idx;
	int64_t lf_pos, dest;

	lf_pos = m->pi[i] + (pos - m->p[i]);
	dest = m->xi[i];

	while (dest + 1 < m->n_runs && m->p[dest + 1] <= lf_pos)
		dest++;
	while (dest > 0 && m->p[dest] > lf_pos)
		dest--;

	*run_idx = dest;
	return lf_pos;
}

int64_t rb3_move_reposition(const rb3_move_t *m, int64_t run_idx, int8_t cc)
{
	return run_idx + m->dist[run_idx * RB3_ASIZE + cc];
}

int64_t rb3_move_step(const rb3_move_t *m, int64_t pos, int64_t *run_idx, int8_t cc)
{
	int64_t i = *run_idx;
	if (m->c[i] != cc) {
		i = rb3_move_reposition(m, i, cc);
		pos = m->p[i];
	}
	*run_idx = i;
	return rb3_move_lf(m, pos, run_idx);
}

/*
 * Binary search for the run containing BWT position pos.
 * Searches the contiguous p[] array (stride 8 bytes, not stride 48).
 */
static int64_t move_find_run(const rb3_move_t *m, int64_t pos)
{
	const int64_t *p = m->p;
	int64_t lo = 0, hi = m->n_runs - 1;
	while (lo < hi) {
		int64_t mid = lo + (hi - lo + 1) / 2;
		if (p[mid] <= pos) lo = mid;
		else hi = mid - 1;
	}
	return lo;
}

int64_t rb3_move_count(const rb3_move_t *m, int len, const uint8_t *pattern)
{
	int64_t lo, hi, *cumrank;
	int i;
	int8_t cc;

	if (len <= 0) return m->bwt_len;
	if (m->n_runs == 0) return 0;

	/* Build cumulative rank table */
	cumrank = RB3_CALLOC(int64_t, ((int64_t)m->n_runs + 1) * RB3_ASIZE);
	{
		int64_t j;
		for (j = 0; j < m->n_runs; j++) {
			int c2;
			for (c2 = 0; c2 < RB3_ASIZE; c2++)
				cumrank[(j + 1) * RB3_ASIZE + c2] = cumrank[j * RB3_ASIZE + c2];
			cumrank[(j + 1) * RB3_ASIZE + m->c[j]] += m->len[j];
		}
	}

	cc = pattern[len - 1];
	if (cc < 0 || cc >= RB3_ASIZE || m->acc[cc] >= m->acc[cc + 1]) {
		free(cumrank);
		return 0;
	}
	lo = m->acc[cc];
	hi = m->acc[cc + 1];

	for (i = len - 2; i >= 0; --i) {
		int64_t lo_run, hi_run, rank_lo, rank_hi;
		cc = pattern[i];
		if (cc < 0 || cc >= RB3_ASIZE || m->acc[cc] >= m->acc[cc + 1]) {
			free(cumrank);
			return 0;
		}

		lo_run = move_find_run(m, lo);
		hi_run = (hi < m->bwt_len) ? move_find_run(m, hi) : m->n_runs - 1;

		rank_lo = cumrank[lo_run * RB3_ASIZE + cc];
		if (m->c[lo_run] == cc)
			rank_lo += lo - m->p[lo_run];

		if (hi >= m->bwt_len) {
			rank_hi = cumrank[m->n_runs * RB3_ASIZE + cc];
		} else {
			rank_hi = cumrank[hi_run * RB3_ASIZE + cc];
			if (m->c[hi_run] == cc)
				rank_hi += hi - m->p[hi_run];
		}

		lo = m->acc[cc] + rank_lo;
		hi = m->acc[cc] + rank_hi;
		if (lo >= hi) {
			free(cumrank);
			return 0;
		}
	}

	free(cumrank);
	return hi - lo;
}

/*
 * .mvi file format — shared header layout (96 bytes):
 *   [0:4]   magic "MVI\1" (v1) or "MVI\2" (v2)
 *   [4:8]   uint32_t flags (reserved)
 *   [8:16]  int64_t n_runs
 *   [16:24] int64_t bwt_len
 *   [24:80] int64_t acc[7]
 *   [80:84] int32_t d
 *   [84:88] uint32_t row_size (v1: 48, v2: 0 = compact)
 *   [88:96] uint64_t checksum
 */

static uint64_t mvi_checksum_v2(const rb3_move_t *m)
{
	int64_t i;
	uint64_t cs = 0;
	for (i = 0; i < m->n_runs; i++) {
		cs ^= (uint64_t)m->xi[i] << 32 | (uint64_t)(uint16_t)m->len[i] << 16 | (uint64_t)(uint8_t)m->c[i];
		cs = (cs << 7) | (cs >> 57); /* rotate */
	}
	return cs;
}

int rb3_move_save(const rb3_move_t *m, const char *fn)
{
	FILE *fp;
	uint32_t flags = 0, row_size = 0; /* row_size=0 signals v2 */
	uint64_t checksum;
	int64_t i;

	fp = fopen(fn, "wb");
	if (fp == 0) return -1;

	/* Check if all lens fit in uint16 (required for v2) */
	{
		int can_v2 = 1;
		for (i = 0; i < m->n_runs; i++) {
			if (m->len[i] > 65535) { can_v2 = 0; break; }
		}
		if (!can_v2) {
			/* Fall back to v1 format for unsplit data */
			rb3_move_row_t row;
			row_size = RB3_MVI_ROW_SIZE_V1;
			fwrite(RB3_MVI_MAGIC_V1, 1, 4, fp);
			fwrite(&flags, 4, 1, fp);
			fwrite(&m->n_runs, 8, 1, fp);
			fwrite(&m->bwt_len, 8, 1, fp);
			fwrite(m->acc, 8, 7, fp);
			fwrite(&m->d, 4, 1, fp);
			fwrite(&row_size, 4, 1, fp);
			/* Compute checksum over v1 rows */
			{
				uint64_t cs = 0;
				for (i = 0; i < m->n_runs; i++) {
					uint64_t *wp;
					memset(&row, 0, sizeof(row));
					row.p = m->p[i]; row.pi = m->pi[i];
					row.xi = m->xi[i]; row.len = m->len[i];
					row.c = m->c[i];
					memcpy(row.dist, &m->dist[i * RB3_ASIZE], sizeof(row.dist));
					wp = (uint64_t *)&row;
					{
						int k;
						for (k = 0; k < (int)(sizeof(row)/sizeof(uint64_t)); k++)
							cs ^= wp[k];
					}
				}
				fwrite(&cs, 8, 1, fp);
			}
			for (i = 0; i < m->n_runs; i++) {
				memset(&row, 0, sizeof(row));
				row.p = m->p[i]; row.pi = m->pi[i];
				row.xi = m->xi[i]; row.len = m->len[i];
				row.c = m->c[i];
				memcpy(row.dist, &m->dist[i * RB3_ASIZE], sizeof(row.dist));
				fwrite(&row, sizeof(row), 1, fp);
			}
			fclose(fp);
			return 0;
		}
	}

	/* V2 compact format */
	checksum = mvi_checksum_v2(m);

	fwrite(RB3_MVI_MAGIC_V2, 1, 4, fp);
	fwrite(&flags, 4, 1, fp);
	fwrite(&m->n_runs, 8, 1, fp);
	fwrite(&m->bwt_len, 8, 1, fp);
	fwrite(m->acc, 8, 7, fp);
	fwrite(&m->d, 4, 1, fp);
	fwrite(&row_size, 4, 1, fp);
	fwrite(&checksum, 8, 1, fp);

	/* Write compact arrays: xi(u32), len(u16), c(i8), dist(i16×6) */
	fwrite(m->xi, sizeof(uint32_t), m->n_runs, fp);
	{
		uint16_t *len16 = RB3_MALLOC(uint16_t, m->n_runs);
		for (i = 0; i < m->n_runs; i++)
			len16[i] = (uint16_t)m->len[i];
		fwrite(len16, sizeof(uint16_t), m->n_runs, fp);
		free(len16);
	}
	fwrite(m->c, sizeof(int8_t), m->n_runs, fp);
	fwrite(m->dist, sizeof(int16_t), m->n_runs * RB3_ASIZE, fp);

	fclose(fp);
	return 0;
}

rb3_move_t *rb3_move_load(const char *fn)
{
	rb3_move_t *m;
	int fd;
	struct stat st;
	uint8_t *base;
	uint32_t flags, row_size;
	uint64_t checksum;
	int is_v2;

	fd = open(fn, O_RDONLY);
	if (fd < 0) return 0;
	if (fstat(fd, &st) != 0 || st.st_size < RB3_MVI_HDR_SIZE) { close(fd); return 0; }

	base = (uint8_t *)mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (base == MAP_FAILED) { close(fd); return 0; }

	/* Check magic */
	is_v2 = (memcmp(base, RB3_MVI_MAGIC_V2, 4) == 0);
	if (!is_v2 && memcmp(base, RB3_MVI_MAGIC_V1, 4) != 0)
		goto fail;

	m = RB3_CALLOC(rb3_move_t, 1);

	/* Parse header */
	memcpy(&flags, base + 4, 4);
	memcpy(&m->n_runs, base + 8, 8);
	memcpy(&m->bwt_len, base + 16, 8);
	memcpy(m->acc, base + 24, 56);
	memcpy(&m->d, base + 80, 4);
	memcpy(&row_size, base + 84, 4);
	memcpy(&checksum, base + 88, 8);

	if (is_v2) {
		/* V2: compact format */
		int64_t i;
		size_t expected = RB3_MVI_HDR_SIZE
			+ (size_t)m->n_runs * sizeof(uint32_t)   /* xi */
			+ (size_t)m->n_runs * sizeof(uint16_t)    /* len */
			+ (size_t)m->n_runs * sizeof(int8_t)      /* c */
			+ (size_t)m->n_runs * RB3_ASIZE * sizeof(int16_t); /* dist */
		if ((size_t)st.st_size != expected) goto fail2;

		move_alloc_arrays(m, m->n_runs);

		/* Read arrays from mmap'd data */
		{
			const uint8_t *ptr = base + RB3_MVI_HDR_SIZE;
			memcpy(m->xi, ptr, m->n_runs * sizeof(uint32_t));
			ptr += m->n_runs * sizeof(uint32_t);

			{
				const uint16_t *len16 = (const uint16_t *)ptr;
				for (i = 0; i < m->n_runs; i++)
					m->len[i] = len16[i];
				ptr += m->n_runs * sizeof(uint16_t);
			}

			memcpy(m->c, ptr, m->n_runs * sizeof(int8_t));
			ptr += m->n_runs * sizeof(int8_t);

			memcpy(m->dist, ptr, m->n_runs * RB3_ASIZE * sizeof(int16_t));
		}

		/* Verify checksum */
		{
			uint64_t computed = mvi_checksum_v2(m);
			if (computed != checksum) goto fail3;
		}

		/* Derive p[] and pi[] */
		move_derive_p(m);
		move_derive_pi(m);

		/* Unmap — we copied the data */
		munmap(base, st.st_size);
		close(fd);

	} else {
		/* V1: legacy 48-byte rows */
		int64_t i;
		const rb3_move_row_t *rows;
		size_t expected;

		if (row_size != RB3_MVI_ROW_SIZE_V1) goto fail2;
		expected = RB3_MVI_HDR_SIZE + (size_t)m->n_runs * row_size;
		if ((size_t)st.st_size != expected) goto fail2;

		rows = (const rb3_move_row_t *)(base + RB3_MVI_HDR_SIZE);

		/* Verify checksum */
		{
			const uint64_t *wp = (const uint64_t *)rows;
			int64_t n_words = m->n_runs * (int64_t)(sizeof(rb3_move_row_t) / sizeof(uint64_t));
			uint64_t cs = 0;
			for (i = 0; i < n_words; i++)
				cs ^= wp[i];
			if (cs != checksum) goto fail2;
		}

		/* Convert v1 rows to separate arrays */
		move_alloc_arrays(m, m->n_runs);
		for (i = 0; i < m->n_runs; i++) {
			m->p[i]   = rows[i].p;
			m->pi[i]  = rows[i].pi;
			m->xi[i]  = (uint32_t)rows[i].xi;
			m->len[i] = rows[i].len;
			m->c[i]   = rows[i].c;
			memcpy(&m->dist[i * RB3_ASIZE], rows[i].dist, RB3_ASIZE * sizeof(int16_t));
		}

		/* Unmap — we copied the data */
		munmap(base, st.st_size);
		close(fd);
	}

	return m;

fail3:
	move_free_arrays(m);
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

int64_t *rb3_move_lcp_thresholds(const rb3_move_t *m, const rb3_lcp_t *lcp)
{
	int64_t *th, i, j;
	if (m == 0 || lcp == 0 || lcp->thresholds == 0) return 0;
	th = RB3_MALLOC(int64_t, m->n_runs);
	j = 0;
	for (i = 0; i < m->n_runs; ++i) {
		while (j + 1 < lcp->n_runs && lcp->run_starts[j + 1] <= m->p[i])
			++j;
		th[i] = lcp->thresholds[j];
	}
	return th;
}

int64_t *rb3_move_lcp_run_map(const rb3_move_t *m, const rb3_lcp_t *lcp)
{
	int64_t *rm, i, j;
	if (m == 0 || lcp == 0) return 0;
	rm = RB3_MALLOC(int64_t, m->n_runs);
	j = 0;
	for (i = 0; i < m->n_runs; ++i) {
		while (j + 1 < lcp->n_runs && lcp->run_starts[j + 1] <= m->p[i])
			++j;
		rm[i] = j;
	}
	return rm;
}

int64_t rb3_move_ms_step(const rb3_move_t *m, const int64_t *run_map,
                         const rb3_lcp_t *lcp,
                         int64_t pos, int64_t *run_idx, int64_t *match_len, int8_t cc)
{
	int64_t i = *run_idx;
	if (cc < 1 || cc >= RB3_ASIZE || m->acc[cc] == m->acc[cc + 1])
		return -1;
	if (m->c[i] == cc) {
		pos = rb3_move_lf(m, pos, run_idx);
		++(*match_len);
	} else {
		int64_t old_lcp_run, new_lcp_run, lo, hi, j, threshold;
		old_lcp_run = run_map[i];
		i = rb3_move_reposition(m, i, cc);
		new_lcp_run = run_map[i];
		if (new_lcp_run < old_lcp_run) {
			if (pos < lcp->tau[old_lcp_run])
				threshold = lcp->lcp_samples[old_lcp_run];
			else
				threshold = lcp->within_min[old_lcp_run];
		} else {
			if (pos >= lcp->tau[old_lcp_run])
				threshold = (old_lcp_run + 1 < lcp->n_runs) ? lcp->lcp_samples[old_lcp_run + 1] : 0;
			else
				threshold = lcp->within_min[old_lcp_run];
		}
		lo = old_lcp_run < new_lcp_run ? old_lcp_run : new_lcp_run;
		hi = old_lcp_run > new_lcp_run ? old_lcp_run : new_lcp_run;
		for (j = lo + 1; j <= hi; ++j)
			if (lcp->lcp_samples[j] < threshold)
				threshold = lcp->lcp_samples[j];
		if (*match_len > threshold)
			*match_len = threshold;
		pos = m->p[i];
		*run_idx = i;
		pos = rb3_move_lf(m, pos, run_idx);
		++(*match_len);
	}
	return pos;
}

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

				{
					int64_t run_idx = move_lcp_find_run(lcp, k);
					int64_t lo_run = run_idx;
					int64_t hi_run = move_lcp_find_run(lcp, l > 0 ? l - 1 : 0);
					while (lo_run > 0 && lcp->lcp_samples[lo_run] >= d)
						--lo_run;
					while (hi_run + 1 < lcp->n_runs && lcp->lcp_samples[hi_run + 1] >= d)
						++hi_run;
					k = lcp->run_starts[lo_run];
					if (hi_run + 1 < lcp->n_runs)
						l = lcp->run_starts[hi_run + 1];
					else
						l = m->bwt_len;
				}

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
 * with SAMPLED cumulative rank table
 ***************************************************/

rb3_bmove_t *rb3_bmove_init(const rb3_move_t *mv)
{
	rb3_bmove_t *bm;
	int64_t i, sample_idx;
	int cc;
	int64_t running[RB3_ASIZE];

	if (mv == 0 || mv->n_runs == 0) return 0;

	bm = RB3_CALLOC(rb3_bmove_t, 1);
	bm->mv = mv;
	bm->n_samples = mv->n_runs / RB3_CUMRANK_SAMPLE + 1;
	bm->cumrank = RB3_CALLOC(int64_t, (bm->n_samples + 1) * RB3_ASIZE);

	/* Build sampled cumrank: store every K-th row's cumulative ranks */
	memset(running, 0, sizeof(running));
	sample_idx = 0;
	/* Sample 0 at row 0: all zeros (already from calloc) */
	for (i = 0; i < mv->n_runs; i++) {
		if (i > 0 && i % RB3_CUMRANK_SAMPLE == 0) {
			sample_idx++;
			for (cc = 0; cc < RB3_ASIZE; cc++)
				bm->cumrank[sample_idx * RB3_ASIZE + cc] = running[cc];
		}
		running[(int)mv->c[i]] += mv->len[i];
	}
	/* Store final totals as last sample (for pos >= bwt_len edge case) */
	sample_idx++;
	for (cc = 0; cc < RB3_ASIZE; cc++)
		bm->cumrank[sample_idx * RB3_ASIZE + cc] = running[cc];
	/* Note: bm->n_samples is the number of sample points (not including the final totals).
	 * Actual entries: n_samples + 1 (the +1 is for the final totals). */

	return bm;
}

void rb3_bmove_destroy(rb3_bmove_t *bm)
{
	if (bm == 0) return;
	free(bm->cumrank);
	free(bm);
}

/*
 * Rank query using sampled cumrank + sequential scan.
 *
 * 1. Binary search on p[] to find run containing pos (stride 8 bytes)
 * 2. Find nearest cumrank sample (run / K)
 * 3. Scan c[] and len[] from sample to run (K sequential accesses)
 * 4. Partial-run adjustment
 */
void rb3_bmove_rank1a(const rb3_bmove_t *bm, int64_t pos, int64_t ok[RB3_ASIZE])
{
	const rb3_move_t *m = bm->mv;
	int64_t run, sample, scan_start, i;
	int cc;

	if (pos <= 0) {
		memset(ok, 0, RB3_ASIZE * sizeof(int64_t));
		return;
	}
	if (pos >= m->bwt_len) {
		/* Use final totals (last entry in sampled cumrank) */
		int64_t last = bm->n_samples;
		for (cc = 0; cc < RB3_ASIZE; cc++)
			ok[cc] = bm->cumrank[last * RB3_ASIZE + cc];
		return;
	}

	/* Binary search on contiguous p[] array */
	run = move_find_run(m, pos);

	/* Find nearest cumrank sample */
	sample = run / RB3_CUMRANK_SAMPLE;
	scan_start = sample * RB3_CUMRANK_SAMPLE;

	/* Start from sampled cumrank */
	for (cc = 0; cc < RB3_ASIZE; cc++)
		ok[cc] = bm->cumrank[sample * RB3_ASIZE + cc];

	/* Scan from sample start to run (sequential, cache-friendly) */
	for (i = scan_start; i < run; i++)
		ok[(int)m->c[i]] += m->len[i];

	/* Partial-run adjustment */
	ok[(int)m->c[run]] += pos - m->p[run];
}

void rb3_bmove_rank2a(const rb3_bmove_t *bm, int64_t k, int64_t l, int64_t ok[RB3_ASIZE], int64_t ol[RB3_ASIZE])
{
	rb3_bmove_rank1a(bm, k, ok);
	rb3_bmove_rank1a(bm, l, ol);
}

void rb3_bmove_extend(const rb3_bmove_t *bm, const rb3_sai_t *ik, rb3_sai_t ok[RB3_ASIZE], int is_back)
{
	int64_t tk[RB3_ASIZE], tl[RB3_ASIZE];
	int cc;

	is_back = !!is_back;
	rb3_bmove_rank1a(bm, ik->x[!is_back], tk);
	rb3_bmove_rank1a(bm, ik->x[!is_back] + ik->size, tl);

	for (cc = 0; cc < RB3_ASIZE; cc++) {
		ok[cc].x[!is_back] = bm->mv->acc[cc] + tk[cc];
		ok[cc].size = (tl[cc] -= tk[cc]);
	}
	ok[0].x[is_back] = ik->x[is_back];
	ok[4].x[is_back] = ok[0].x[is_back] + tl[0];
	ok[3].x[is_back] = ok[4].x[is_back] + tl[4];
	ok[2].x[is_back] = ok[3].x[is_back] + tl[3];
	ok[1].x[is_back] = ok[2].x[is_back] + tl[2];
	ok[5].x[is_back] = ok[1].x[is_back] + tl[1];
}

static inline void bmove_set_intv(const rb3_bmove_t *bm, int cc, rb3_sai_t *ik)
{
	const int64_t *acc = bm->mv->acc;
	ik->x[0] = acc[cc];
	ik->size = acc[cc + 1] - acc[cc];
	ik->x[1] = acc[rb3_comp(cc)];
	ik->info = 0;
}

static void bmove_sai_reverse(rb3_sai_t *a, int64_t l)
{
	int64_t i;
	rb3_sai_t t;
	for (i = 0; i < l>>1; ++i)
		t = a[i], a[i] = a[l - 1 - i], a[l - 1 - i] = t;
}

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
	for (i = x + 1, curr->n = 0; i < len; ++i) {
		int cc = rb3_comp(q[i]);
		rb3_bmove_extend(bm, &ik, ok, 0);
		if (ok[cc].size != ik.size) {
			Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
			curr->a[curr->n++] = ik;
			if (ok[cc].size < min_occ) break;
		}
		ik = ok[cc]; ik.info = i + 1;
	}
	if (i == len) {
		Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
		curr->a[curr->n++] = ik;
	}
	bmove_sai_reverse(curr->a, curr->n);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) {
		int cc = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < (int64_t)prev->n; ++j) {
			rb3_sai_t *p = &prev->a[j];
			rb3_bmove_extend(bm, p, ok, 1);
			if (cc == 0 || ok[cc].size < min_occ) {
				if (curr->n == 0 && (int32_t)p->info - i - 1 >= min_len && (mem->n == oldn || i + 1 < mem->a[mem->n-1].info>>32)) {
					rb3_sai_t *q2;
					Kgrow(km, rb3_sai_t, mem->a, mem->n, mem->m);
					q2 = &mem->a[mem->n++];
					*q2 = *p; q2->info |= (int64_t)(i + 1)<<32;
				}
			} else if (curr->n == 0 || ok[cc].size != curr->a[curr->n-1].size) {
				ok[cc].info = p->info;
				Kgrow(km, rb3_sai_t, curr->a, curr->n, curr->m);
				curr->a[curr->n++] = ok[cc];
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

static int64_t rb3_bmove_smem1_TG(void *km, const rb3_bmove_t *bm, int64_t min_occ, int64_t min_len, int64_t len, const uint8_t *q, int64_t x, rb3_sai_v *mem, int32_t check_long)
{
	int64_t i, j;
	rb3_sai_t ik, ok[6];

	assert(len <= INT32_MAX);
	if (len - x < min_len) return len;
	bmove_set_intv(bm, q[x + min_len - 1], &ik);
	for (i = x + min_len - 2; i >= x; --i) {
		int cc = q[i];
		rb3_bmove_extend(bm, &ik, ok, 1);
		if (ok[cc].size < min_occ) break;
		ik = ok[cc];
	}
	if (i >= x) return i + 1;
	if (check_long) return -1;
	for (j = x + min_len; j < len; ++j) {
		int cc = rb3_comp(q[j]);
		rb3_bmove_extend(bm, &ik, ok, 0);
		if (ok[cc].size < min_occ) break;
		ik = ok[cc];
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
	for (i = j - 1; i > x; --i) {
		int cc = q[i];
		rb3_bmove_extend(bm, &ik, ok, 1);
		if (ok[cc].size < min_occ) break;
		ik = ok[cc];
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
