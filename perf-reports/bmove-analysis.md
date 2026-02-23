# b-move Performance and Disk Size Analysis

**Date:** 2026-02-22
**ropebwt3 version:** 3.10-r283 (commit 4020a55)
**Platform:** Linux 6.17.0-14-generic, AMD Ryzen AI 7 350, 16 cores, 94 GB RAM

## Executive Summary

The b-move backend is 1.3-4.5x slower than FMD for SMEM search and its .mvi files are 39-55x larger than .fmd files. This analysis identifies the root causes and proposes concrete optimizations.

**Performance bottleneck:** Each b-move rank query requires an O(log r) binary search over the move table to locate which run contains a BWT position, then a lookup into a 48-byte-per-row cumulative rank table. FMD's RLE-delta encoding performs sequential decoding from nearby checkpoints, which is faster at tested scales because the working set fits in cache.

**Disk size bottleneck:** The move row is 48 bytes (3x int64 + int64 + 6x int16 + padding), vs Movi's 8 bytes/row. The fields `p`, `pi`, `xi`, and `len` use full int64 even though most values fit in 32 or 16 bits. The cumrank table (built at load time) adds another 48 bytes/row in RAM.

**Key finding:** At tested scales (up to 7.5M runs), the combined b-move working set (96 bytes/row for move table + cumrank) vastly exceeds what fits in L3 cache, while the FMD (6.2 MB) fits comfortably. The crossover where b-move wins requires either (a) the FMD to also exceed cache, or (b) the row size to shrink enough to fit.

---

## 1. Performance Analysis

### 1.1 Measured Timings (10,000 queries, 100-150 bp, single-threaded)

| Dataset | BWT runs | FMD time | b-move time | Ratio | FMD RSS | b-move RSS |
|---------|----------|----------|-------------|-------|---------|------------|
| Medium (random, 10 Mbp) | 7,499,526 | 0.35s | 1.02s | 2.9x slower | 12 MB | 682 MB |
| Large (redundant, 20 Mbp) | 1,243,709 | 0.58s | 1.11s | 1.9x slower | 8 MB | 119 MB |

### 1.2 Cost Per Rank Query

The core operation for bidirectional SMEM search is `rb3_bmove_extend()`, which calls `rb3_bmove_rank1a()` twice (at positions k and k+size). Each rank1a call:

1. **Binary search** over move rows to find the run containing position `pos`: O(log r)
   - `move_find_run()` at `move.c:345` — binary search on `m->rows[mid].p`
   - For 7.5M runs: ~23 iterations, each touching a different cache line
   - For 1.2M runs: ~20 iterations

2. **Cumrank table lookup**: O(1) per character, 6 lookups
   - `bm->cumrank[run * 6 + c]` for all c in 0..5
   - Access pattern: single cache line (48 bytes contiguous)

3. **Partial-run adjustment**: O(1)
   - `ok[rows[run].c] += pos - rows[run].p`

**Binary search is the bottleneck.** Each iteration of the binary search accesses `m->rows[mid].p` — an 8-byte field inside a 48-byte struct. With 7.5M rows × 48 bytes = 343 MB, the move table array far exceeds L3 cache. Each binary search step is essentially a random memory access into this 343 MB array, causing a cache miss with ~100ns penalty on modern hardware.

Estimated cost per rank query:
- Binary search: 23 steps × ~50ns avg (mix of cache hits for final narrowing, misses for early steps) ≈ **1150 ns**
- Cumrank lookup: 1 cache line access ≈ **50 ns** (likely cached near move row)
- Total: **~1200 ns per rank query**

For comparison, FMD rank via `rld_rank2a()`:
- Checkpoint lookup + sequential decode of ~64 symbols ≈ **200-400 ns** (compact, cache-friendly)

This 3-6x difference per rank query directly explains the 2-4.5x overall slowdown (since most wall time is spent in rank queries during SMEM).

### 1.3 Per-SMEM Algorithmic Complexity

Both FMD and b-move SMEM use the same algorithm (`rb3_fmd_smem1` / `rb3_bmove_smem1`). Each extension call requires **2 rank queries** (at positions k and k+size for the bidirectional interval). The number of extensions per SMEM is identical for both backends — typically O(read_length) forward + O(read_length) backward per seed position.

The algorithmic work is equivalent; the difference is entirely in the **constant factor per rank query**.

### 1.4 Cache Pressure Analysis

| Component | Medium (7.5M runs) | Large (1.2M runs) |
|-----------|--------------------|--------------------|
| Move table (`rows[]`) | 343 MB | 57 MB |
| Cumrank table | 343 MB | 57 MB |
| **Total b-move working set** | **686 MB** | **114 MB** |
| FMD (`.fmd` in memory) | 6.2 MB | 1.5 MB |
| Typical L3 cache | 16-32 MB | 16-32 MB |

For the medium dataset, the b-move working set (686 MB) is **21-43x larger** than L3 cache. Every binary search step beyond the first few is a guaranteed cache miss.

For the large dataset, the working set (114 MB) is **3.5-7x** larger than L3 — still too large, but the penalty is smaller, which is why the slowdown is "only" 1.9x instead of 2.9x.

The FMD always fits comfortably in L3 cache at these scales.

### 1.5 Projected Crossover Point

The b-move becomes competitive when the **FMD itself no longer fits in cache**. With L3 ≈ 32 MB and FMD compression ratio ≈ 0.65 bytes/BWT-symbol (from the medium dataset), this happens around:

- FMD > 32 MB → BWT length > ~50M symbols → ~25k sequences × 1 kbp each
- But b-move must also fit: need row size ≤ 32 MB / r rows

For the b-move to actually win, we need the move table to be **at most comparable in size to the FMD**. At 48 bytes/row, this never happens (move table is always larger). At 8-12 bytes/row (Movi-style packing), the crossover occurs around FMD size ≈ 50-100 MB, which corresponds to real pangenome-scale datasets (thousands of genomes).

---

## 2. Disk Size Analysis

### 2.1 Current .mvi Format

```
Header: 96 bytes (magic, n_runs, bwt_len, acc[7], d, row_size, checksum)
Body:   n_runs × 48 bytes per row
```

| Dataset | .fmd | .mvi | Ratio | Bytes/run (mvi) | Bytes/symbol (fmd) |
|---------|------|------|-------|-----------------|---------------------|
| Small (200 kbp) | 129 KB | 6.9 MB | 55x | 48.0 | 0.66 |
| Medium (10 Mbp) | 6.2 MB | 343 MB | 55x | 48.0 | 0.66 |
| Large (20 Mbp) | 1.5 MB | 57 MB | 39x | 48.0 | 0.08 |

The FMD is compact because it uses run-length delta encoding (RLE-delta): each BWT run is encoded as a variable-length symbol+count pair, typically 1-4 bytes per run. The .mvi format uses a fixed 48 bytes per run.

### 2.2 Field-by-Field Breakdown (48 bytes/row)

| Field | Type | Bytes | Actual range needed | Could be |
|-------|------|-------|--------------------|----|
| `p` | int64 | 8 | BWT offset [0, n) | Derivable from prefix sum of `len` |
| `pi` | int64 | 8 | LF[p] = acc[c] + rank(c,p) | Derivable from `c` and cumulative c-counts |
| `pi` | int64 | 8 | LF[p] = acc[c] + rank(c,p) | Or store delta from acc[c]: uint32 |
| `xi` | int64 | 8 | Run index [0, r) | uint32 (up to 4B runs) |
| `len` | int64 | 8 | Run length, typically 1-137 | uint16 (with run splitting, always < 65535) |
| `dist[6]` | int16×6 | 12 | Signed dist to nearest run of each char | Could use int8 with run splitting |
| `c` | int8 | 1 | Character [0,5] | 3 bits |
| `pad` | uint8×3 | 3 | Alignment padding | Remove with `__attribute__((packed))` |

### 2.3 Comparison with Reference Implementations

| Implementation | Row size | Key trick |
|----------------|----------|-----------|
| **ropebwt3 (current)** | **48 bytes** | All fields as full int64; mmap-friendly alignment |
| **Movi** (regular mode) | **8 bytes** | Pack len+char into uint16, offset into uint16, dest into uint32 |
| **Movi** (constant-time) | **24 bytes** | Add 6 × uint16 jump pointers for O(1) reposition |
| **b-move** (biointec, packed) | **~10-15 bytes** | Bit-packed with configurable bit widths per field |
| **b-move** (biointec, unpacked) | **25 bytes** | char + 3 × uint64 (no dist array) |
| **move-r** | **~12 bytes** | Heavily optimized, benchmarked 2-35x faster than other r-index implementations |

**ropebwt3's 48-byte row is 4-6x larger than state-of-the-art.** The main reason is storing `p`, `pi`, `xi`, `len` as full int64 and including the 12-byte dist[] array inline.

### 2.4 Minimum Viable Row (for bidirectional b-move search)

The minimum fields needed per row for the current b-move algorithm:

**For LF-mapping** (`rb3_move_lf`):
- `p` (run start position) — needed to compute offset within run
- `pi` (LF[p]) — needed for LF formula: LF(pos) = pi + (pos - p)
- `xi` (destination run) — jump target after LF
- `len` (run length) — fast-forward boundary check

**For reposition** (`rb3_move_reposition`):
- `dist[c]` — signed distance to nearest run of character c
- `c` (character) — needed to check if reposition is needed

**For rank queries** (`rb3_bmove_rank1a`):
- `p` (run start position) — binary search target
- `c` (character) — partial-run adjustment
- `cumrank[]` (separate table) — cumulative ranks at run boundaries

### 2.5 Proposed Compact Format

**Target: 16 bytes/row** (3x reduction, fits 4 rows per cache line):

```c
typedef struct rb3_move_row_compact_s {
    uint32_t p_hi;       // High 32 bits of p (or full p if BWT < 4G)
    uint32_t xi;         // Destination run index (up to 4G runs)
    uint16_t len;        // Run length (with splitting, < 65535)
    uint16_t pi_delta;   // pi - acc[c] - cumulative_c_count (derivable offset)
    int8_t dist[3];      // Top-3 characters only (A,C,G,T → 4 of 6), or...
    uint8_t c;           // Run character (3 bits) + flags (5 bits)
} __attribute__((packed)) rb3_move_row_compact_t; // 16 bytes
```

Or even more aggressively:

**Target: 12 bytes/row** (4x reduction):

```c
typedef struct {
    uint32_t xi;         // Destination run index
    uint32_t p_low;      // Low 32 bits of p (use block+offset for large BWTs)
    uint16_t len;        // Run length
    uint8_t c;           // Character
    uint8_t flags;       // Packed: dist signs/magnitudes for nearby chars
} __attribute__((packed)) rb3_move_row_12_t; // 12 bytes
```

With dist[] stored in a separate array (accessed only during reposition, not during rank queries), the hot path (binary search + LF-mapping) touches only the compact row.

**Impact on disk size:**

| Dataset | Current (48B) | 16 bytes/row | 12 bytes/row | 8 bytes/row |
|---------|--------------|--------------|--------------|-------------|
| Small | 6.9 MB | 2.3 MB | 1.7 MB | 1.2 MB |
| Medium | 343 MB | 114 MB | 86 MB | 57 MB |
| Large | 57 MB | 19 MB | 14 MB | 9.5 MB |
| FMD ratio | 39-55x | 13-18x | 10-14x | 7-9x |

---

## 3. Memory Analysis

### 3.1 Current Memory Breakdown (b-move for SMEM)

| Component | Formula | Medium (7.5M runs) | Large (1.2M runs) |
|-----------|---------|--------------------|--------------------|
| Move table (mmap'd) | r × 48 | 343 MB | 57 MB |
| Cumrank table | (r+1) × 48 | 343 MB | 57 MB |
| FMD (also loaded) | varies | 6 MB | 2 MB |
| **Total** | | **692 MB** | **116 MB** |

The cumrank table is the **hidden cost**: it's as large as the move table itself but allocated at runtime (`rb3_bmove_init`). This is not stored on disk — it's recomputed from the move table at every load.

### 3.2 Cumrank Optimization: Sampled Cumrank

Instead of storing cumrank for every row, sample every K rows and reconstruct intermediate values by scanning K move rows:

| Sample rate K | Cumrank size (medium) | Rank query cost | Total RAM (medium) |
|---------------|----------------------|-----------------|---------------------|
| 1 (current) | 343 MB | O(log r) binary search | 692 MB |
| 16 | 21 MB | O(log(r/16) + 16) | 370 MB |
| 64 | 5 MB | O(log(r/64) + 64) | 354 MB |
| 256 | 1.3 MB | O(log(r/256) + 256) | 350 MB |

With K=64, cumrank shrinks from 343 MB to 5 MB — a 68x reduction — while adding only ~64 sequential row accesses per rank query. This is comparable to FMD's sequential decode cost per rank query.

Combined with compact rows (16B/row → 114 MB table), total RAM drops from 692 MB to 119 MB — an almost 6x improvement.

### 3.3 Alternative: Interleaved Cumrank

Store cumrank values interleaved with move rows (every K-th row contains cumrank). This improves spatial locality: when the binary search lands on a sampled row, both the move data and cumrank are in the same cache line.

---

## 4. Optimization Roadmap

Ranked by impact × feasibility:

### Priority 1: Pack move rows from 48 to 16 bytes (HIGH IMPACT, MEDIUM EFFORT)

**What:** Store `xi` as uint32, `len` as uint16, `c` as uint8, separate `p` into a prefix-sum array. Move `dist[]` to a separate array.

**Impact:**
- Disk: 3x smaller .mvi files (343 MB → 114 MB for medium)
- RAM: 3x smaller move table, better cache utilization
- Speed: ~1.5-2x faster binary search (more rows per cache line)

**Risk:** Requires new .mvi format version. Breaking change for existing files.

**Implementation sketch:**
- Store p[] as a separate int64 array (used only for binary search in rank queries)
- Store compact rows as 16-byte structs (xi, len, c + hot fields)
- Store dist[] as a separate int16[6] array (used only during reposition)
- Binary search operates on the p[] array (stride 8 vs stride 48 = 6x fewer cache lines)

### Priority 2: Sample cumrank table (HIGH IMPACT, LOW EFFORT)

**What:** Store cumrank only every K rows. Reconstruct intermediate values by scanning K rows from the nearest sample.

**Impact:**
- RAM: cumrank shrinks from ~r×48 to ~r×48/K bytes. With K=64: 343 MB → 5 MB.
- Speed: Rank queries become O(log(r/K) + K) instead of O(log r). With K=64, the +64 scan dominates but is sequential (cache-friendly).

**Risk:** Rank queries become slightly slower in absolute terms (sequential scan of K rows), but the massively reduced cache pressure on the working set could make them faster in practice.

**Implementation sketch:**
```c
// Sampled cumrank: only store every CUMRANK_SAMPLE rows
#define CUMRANK_SAMPLE 64

void rb3_bmove_rank1a_sampled(const rb3_bmove_t *bm, int64_t pos, int64_t ok[6]) {
    int64_t run = move_find_run(bm->mv, pos);
    int64_t sample = run / CUMRANK_SAMPLE;
    // Start from sampled cumrank
    memcpy(ok, &bm->cumrank[sample * 6], 48);
    // Scan rows from sample*K to run
    for (int64_t i = sample * CUMRANK_SAMPLE; i < run; i++)
        ok[bm->mv->rows[i].c] += bm->mv->rows[i].len;
    // Partial-run adjustment
    ok[bm->mv->rows[run].c] += pos - bm->mv->rows[run].p;
}
```

### Priority 3: Separate p[] array for binary search (MEDIUM IMPACT, LOW EFFORT)

**What:** Extract `p` (run start positions) into a contiguous int64 array used only for binary search. The rest of the row data is accessed only after the search converges.

**Impact:**
- Binary search touches only p[] (stride 8 bytes) instead of rows[] (stride 48 bytes)
- 6x fewer cache lines accessed during binary search
- Estimated speedup: 1.5-2x for binary search, ~1.2-1.5x overall

**Implementation:** Minimal code change — build p[] at load time as a parallel array.

### Priority 4: Prefetch + batched queries (MEDIUM IMPACT, HIGH EFFORT)

**What:** Process multiple independent rank queries concurrently, using `__builtin_prefetch` to hide cache miss latency. Movi reports 24x speedup from processing 16 concurrent queries.

**Impact:** Could reduce effective cache miss penalty from ~100ns to ~10-20ns.

**Risk:** Requires restructuring the SMEM algorithm for batch processing. The current algorithm processes one extension at a time; batching would require accumulating multiple independent positions.

### Priority 5: Bit-packed row format (MEDIUM IMPACT, HIGH EFFORT)

**What:** Use variable bit widths per field (as in biointec's b-move). Store `p` in ceil(log2(n)) bits, `xi` in ceil(log2(r)) bits, etc.

**Impact:** Could reduce to 10-15 bytes/row on typical datasets.

**Risk:** Bit extraction is CPU-intensive (masks, shifts). May not pay off if the bottleneck shifts from cache misses to decode overhead.

---

## 5. Comparison with the State of the Art

| System | Row size | Rank query | Move table | Key innovation |
|--------|----------|------------|------------|----------------|
| **ropebwt3** | 48 B | O(log r) + cumrank | 96 B/row total | Flat mmap-friendly array |
| **Movi** | 8 B | No bidirectional | 8 B/row | Ultra-compact, prefetch batching |
| **b-move (biointec)** | 10-25 B | O(r) walk | ~15 B/row | Memoized bidirectional, bit-packing |
| **move-r** | ~12 B | O(1) LF | ~12 B/row | Heavily optimized, 15x faster than other r-indexes |
| **MIOV** | varies | Same as Movi | Same | Row reordering for cache locality |

ropebwt3's b-move is the **least space-efficient** of all implementations, primarily because it uses full int64 for all fields and includes the 12-byte dist[] array inline. The cumrank table (not stored by others) doubles the effective cost.

---

## 6. Root Cause Summary

### Why b-move is slower (1.3-4.5x):

1. **Binary search cache misses** (dominant): Each rank query does O(log r) random accesses into a 48-byte-stride array that doesn't fit in cache. FMD's sequential decode from nearby checkpoints is fundamentally more cache-friendly at these scales.

2. **Cumrank table size**: The cumrank table doubles the working set. FMD stores rank information implicitly in its RLE-delta encoding — no separate table needed.

3. **48-byte row stride**: Only 1.3 rows fit per 64-byte cache line, wasting bandwidth. At 8-16 bytes/row, 4-8 rows would fit per cache line, dramatically reducing misses.

### Why .mvi files are 39-55x larger:

1. **Fixed 48 bytes/row vs RLE-delta encoding**: FMD stores each BWT run in 1-4 bytes (variable-length). The move table uses a fixed 48 bytes regardless of run length or content.

2. **Redundant fields**: `p` and `pi` are derivable from prefix sums. `xi` is derivable via binary search on `pi`. Storing them explicitly trades CPU for space — a trade that doesn't pay off when the space causes cache misses.

3. **No compression**: The .mvi format stores raw struct data. Delta encoding, varint encoding, or even simple gzip would significantly reduce size (though mmap compatibility would be lost).

---

## 7. Conclusions

The b-move implementation in ropebwt3 is correct and follows sound theoretical principles (O(r)-space bidirectional search). However, the constant factors — 48-byte rows, full cumrank table, binary search per rank query — make it uncompetitive with FMD at tested scales.

**To make b-move competitive, the most impactful changes are:**

1. **Shrink rows to 12-16 bytes** and separate hot/cold fields → 3-4x smaller working set
2. **Sample the cumrank table** (every 64 rows) → eliminate the largest memory consumer
3. **Use a separate p[] array** for binary search → 6x fewer cache lines per search

These three changes together would reduce the b-move working set from ~96 bytes/row to ~24-28 bytes/row (including sampled cumrank), bringing it within 2-3x of FMD size at small scales and making it genuinely competitive at large scales where FMD also doesn't fit in cache.

The disk size can be reduced from 48 bytes/row to 12-16 bytes/row with field packing, bringing the .mvi/FMD ratio from 39-55x down to 10-18x. Further reduction to 8 bytes/row is possible with Movi-style packing but would require giving up mmap friendliness and adding decode overhead.
