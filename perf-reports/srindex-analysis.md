# SR-Index Deep Analysis: Disk Size and Performance

**Date:** 2026-02-22
**Commit:** 4020a55 (base), worktree `.claude/worktrees/srindex-study`
**Platform:** Linux 6.17.0-14-generic

## Executive Summary

The SR-index (.sri) is 80-200x larger than the sampled suffix array (.ssa) because it stores
six 64-bit arrays at O(r) or O(n/s) scale using raw int64_t, while the SSA packs one uint64_t
per sample with a sparse sampling rate. The two biggest problems are:

1. **For s=1, sub_pos/sub_sa are 100% identical copies of run_pos/run_sa** (wastes 33% of the file)
2. **All arrays use 64-bit integers regardless of actual value range** (18-25 bits suffice for tested datasets)

Two optimizations were implemented and tested:
- Eliminating the s=1 sub array duplication → **33% smaller s=1 index** (59.7→39.8 MB)
- Adding a bitvector for O(1) sub_pos lookup in `locate_one` → faster toehold resolution

## 1. SRI File Format: What's Stored

The SRI file stores these arrays as raw 64-bit integers:

```
Header (48 bytes):
  magic[4] + s(4) + m(8) + n(8) + n_runs(8) + n_samples(8) + n_sub(8)

Arrays:
  phi_sa[n_runs]      - sorted SA values at BWT run starts (phi breakpoints)
  phi_da[n_runs]      - phi values at breakpoints: phi(phi_sa[i]) = phi_da[i]
  run_pos[n_samples]  - BWT position of last char in each run (toehold support)
  run_sa[n_samples]   - SA value at each run_pos (toehold values)
  sub_pos[n_sub]      - BWT positions for subsampled SA (locate_one targets)
  sub_sa[n_sub]       - SA values at sub_pos positions
  cum_len[m+1]        - cumulative sequence lengths (multi-string support)
  text_order_sid[m]   - sentinel BWT positions in text order
```

Where: n_samples = n_runs (one sample per run boundary end), n_sub = n_runs for s=1 or ~n/s for s>1.

## 2. Size Breakdown

### Dataset: large (repetitive, r/n=0.062)
n=20,000,400, r=1,243,709, m=400

| Component | s=1 (bytes) | s=1 (%) | s=8 (bytes) | s=8 (%) |
|-----------|-------------|---------|-------------|---------|
| phi_sa+phi_da | 19,899,344 | 33.3% | 19,899,344 | 24.9% |
| run_pos+run_sa | 19,899,344 | 33.3% | 19,899,344 | 24.9% |
| sub_pos+sub_sa | 19,899,344 | 33.3% | 40,006,400 | 50.1% |
| cum_len+tosid | 6,408 | 0.0% | 6,408 | 0.0% |
| **Total** | **59,704,488** | | **79,811,544** | |

**Key observation for s=1**: sub_pos and sub_sa are byte-for-byte identical to run_pos and run_sa.
This is because `build_subsampled()` does `memcpy(sub_pos, run_pos, ...)` for s≤1.

### Dataset: small (non-repetitive, r/n=0.752)
n=200,200, r=150,569, m=200

| Index | s | Size (bytes) | Size (human) | SSA comparison |
|-------|---|-------------|-------------|----------------|
| SSA | 8 | 7,884 | 7.7 KB | 1.0x |
| SRI | 1 | 7,230,568 | 6.9 MB | 917x |
| SRI | 8 | 5,224,664 | 5.0 MB | 663x |
| SRI | 64 | 4,872,664 | 4.6 MB | 618x |

For non-repetitive data (r≈0.75n), the SRI is 600-900x larger than SSA.

## 3. Array Compression Analysis

Diagnostic analysis of each array's value distribution:

### phi_sa (sorted, delta-encodable)

| Dataset | n_runs | max_value | bits_needed | avg_delta | max_delta | est. delta size |
|---------|--------|-----------|-------------|-----------|-----------|-----------------|
| small (r/n=0.75) | 150,569 | 200,199 | 18 | 1.3 | 12 | 41.6 KB vs 1,176 KB raw |
| large (r/n=0.06) | 1,243,709 | 20,000,399 | 25 | 16.1 | 864 | 500 KB vs 9,716 KB raw |

phi_sa is sorted → deltas are small positive integers. Delta encoding achieves **19-28x compression**
over raw 64-bit storage.

### run_pos (sorted, delta-encodable)

Similar to phi_sa: sorted BWT positions with small deltas (= BWT run lengths).
Delta encoding achieves **19-28x compression**.

### phi_da, run_sa, sub_sa (unsorted, not directly compressible)

These contain arbitrary SA values. Not amenable to delta encoding without a rearrangement.
However, all values fit in 18-25 bits, so **bit-width reduction** to 32-bit would halve them.

### 32-bit feasibility

For n < 2^32 (≈4.3 billion), ALL values in ALL arrays fit in 32 bits:

| Dataset | n | max_value | Current (64-bit) | 32-bit | Savings |
|---------|---|-----------|-------------------|--------|---------|
| small | 200K | 200,199 | 7.2 MB | 3.6 MB | 50% |
| large | 20M | 20,000,399 | 79.8 MB | 39.9 MB | 50% |

## 4. Comparison with Reference SR-Index

The reference implementation (Cobas et al., [duscob/sr-index](https://github.com/duscob/sr-index))
uses SDSL (Succinct Data Structure Library) for compressed storage:

- **Phi breakpoints**: stored in an `sd_vector<>` (sparse bitvector with Elias-Fano encoding)
  with rank/select support, not as raw sorted arrays
- **SA samples**: stored in `sdsl::int_vector<>` with automatically computed bit-width
  (uses ⌈log₂ r⌉ bits per entry, not 64)
- **Toehold**: computed on-demand using BWT rank/select operations

Our implementation uses **no compression at all** - every value is a raw int64_t. This explains
the dramatic size difference. The SDSL approach achieves near-theoretical space bounds while
our approach trades space for simplicity and query speed (no decompression overhead).

## 5. Performance Analysis

### Why is s=1 slower than s=4/s=8?

`locate_one` walks LF from a BWT position until it finds a position in sub_pos[].

- **For s=1**: sub_pos = run boundary endpoints, covering r out of n BWT positions.
  Probability of hitting one ≈ r/n. Expected walk: **n/r steps** (≈16 for large dataset).
- **For s=k (k>1)**: sub_pos contains positions where SA % k == 0.
  Since SA decreases by 1 per LF step, max walk is k-1, expected **(k-1)/2 steps**.

| s value | Expected locate_one walk | Measured (us/call) |
|---------|-------------------------|--------------------|
| 1 | n/r ≈ 16 | 9.6 |
| 4 | 1.5 | ~1.0 (est.) |
| 8 | 3.5 | 0.9 |
| 16 | 7.5 | ~3.5 (est.) |
| 64 | 31.5 | ~15 (est.) |

s=1 is **10.5x slower** than s=8 for locate_one because it takes ~16 LF steps per lookup
vs ~4 for s=8. Despite "no subsampling", s=1 has fewer sample points to find than s=4 or s=8.

### locate_all: phi chain performance

For locate_all, the phi chain enumeration is the same for all s values (O(log r) binary search
per phi step). The only difference is toehold resolution (first locate_one call). For large
intervals, toehold resolution is amortized and s values converge:

| Metric | s=1 | s=8 |
|--------|-----|-----|
| locate_one | 9.6 us | 0.9 us |
| locate_all (avg 25 occ) | 17.3 us | 6.8 us |
| Per-occurrence after toehold | ~0.3 us | ~0.2 us |

The per-occurrence cost after toehold resolution is similar (~0.2-0.3 us) because it's
dominated by the O(log r) phi binary search, which is identical for all s values.

### Why SRI is faster than SSA for high-p

SSA uses the BWT interval splitting approach: it recursively splits the [lo,hi) interval
by LF-mapping until it hits sampled positions. Each level does O(ASIZE) rank queries per
split. For large intervals, this creates O(interval * n/2^ss) total work.

SRI uses phi chain: one locate_one call + (interval-1) phi evaluations. Each phi evaluation
is a single O(log r) binary search, independent of n. For highly repetitive data (r << n),
this is dramatically faster.

## 6. Implemented Optimizations

### 6a. Eliminate s=1 sub array duplication (IMPLEMENTED)

**Problem**: For s=1, `build_subsampled()` allocates new arrays and memcpy's run_pos/run_sa
into sub_pos/sub_sa. This wastes 33% of disk space and memory.

**Fix**: For s=1, sub_pos/sub_sa now alias run_pos/run_sa (shared pointers). The serialization
format is bumped to v2 (magic "SRI\2") which writes n_sub=0 for s=1 to indicate the alias.
The restore function detects this and reconstructs the alias. Backward compatibility with v1
files is maintained.

**Impact**:
- s=1 disk: 59.7 MB → 39.8 MB (**33% smaller**)
- s=1 memory: -20 MB (no redundant allocation)
- s>1: no change

### 6b. Bitvector for sub_pos lookup in locate_one (IMPLEMENTED)

**Problem**: `locate_one` does O(log n_sub) binary search at every LF step to check if the
current position is in sub_pos[]. For s=8 with n_sub=2.5M, that's 21 comparisons per step.

**Fix**: Build a bitvector of n bits (n/8 bytes) marking all sub_pos positions. In locate_one,
check the bitvector first (O(1) single-word AND), and only do the binary search on a hit.

**Impact**:
- Bitvector size: n/8 bytes (2.5 MB for n=20M, vs 20 MB for sub_pos at s=8)
- Each non-hit LF step: 1 memory access (bitvector word) instead of ~21 (binary search)
- Amortized improvement depends on hit rate; most significant for large s where hits are rare

## 7. Optimization Opportunities NOT Yet Implemented

### 7a. 32-bit integer mode (HIGH IMPACT: ~50% size reduction)

For n < 2^32, switch all arrays from int64_t to int32_t. This would halve the index size:
- s=8 large: 79.8 MB → 39.9 MB
- s=1 large (with dedup): 39.8 MB → 19.9 MB

**Complexity**: Moderate. Requires either templating the code or maintaining parallel 32/64-bit
paths. The file format would need a flag byte to indicate width.

### 7b. Delta-encode phi_sa and run_pos (HIGH IMPACT on sorted arrays)

phi_sa and run_pos are sorted arrays with small deltas:
- phi_sa deltas: avg 1.3 (non-repetitive) to 16.1 (repetitive), max 12-864
- run_pos deltas: avg 1.3 to 16.1, max 10-550

Using variable-length delta encoding (e.g., Elias-gamma or byte-aligned vbyte):
- phi_sa: 9.7 MB → 0.5 MB (**19x compression**)
- run_pos: 9.7 MB → 0.5 MB (**19x compression**)

**Trade-off**: Binary search on compressed arrays requires an auxiliary index (e.g., sampling
every 64th prefix sum), adding O(r/64) overhead but providing O(log r + 64) lookup. The
reference sr-index uses SDSL's sd_vector with Elias-Fano for this exact purpose.

**Complexity**: High. Requires implementing or importing compressed vector with rank support.

### 7c. Eliminate run_pos redundancy

The run_pos array stores BWT positions of run ENDS. These are determined by run boundaries,
which are also implicitly captured in phi_sa (sorted SA values at run STARTS). In principle,
run_pos could be reconstructed from the BWT structure, but it would require an additional
scan. In practice, maintaining run_pos for O(log r) toehold lookup is worthwhile.

### 7d. Bit-pack phi_da and run_sa to min_bits

phi_da and run_sa contain arbitrary SA values in [0, n). With n < 2^25 (large dataset),
they need only 25 bits per entry instead of 64. A fixed-width bit-packed array at
⌈log₂ n⌉ bits per entry would save:
- phi_da: 9.9 MB → 3.9 MB (25/64 ratio)
- run_sa: 9.9 MB → 3.9 MB

**Complexity**: Low-moderate. Fixed-width bit-packing is simpler than variable-length encoding
but requires careful bit manipulation.

## 8. Size Projection After All Optimizations

For the large dataset (n=20M, r=1.24M, s=8):

| Optimization | Current | After | Savings |
|-------------|---------|-------|---------|
| Baseline | 79.8 MB | - | - |
| + s=1 dedup (only for s=1) | - | - | (s=8 unchanged) |
| + 32-bit mode | 79.8 MB | 39.9 MB | 50% |
| + delta-encode phi_sa+run_pos | 39.9 MB | 29.9 MB | 25% more |
| + bit-pack phi_da+run_sa (25 bits) | 29.9 MB | 22.2 MB | 26% more |
| + bit-pack sub arrays (25 bits) | 22.2 MB | 14.4 MB | 35% more |
| **All combined** | **79.8 MB** | **~14 MB** | **82%** |

For s=1 (after dedup):

| Optimization | Current | After | Savings |
|-------------|---------|-------|---------|
| Baseline (with dedup) | 39.8 MB | - | - |
| + 32-bit mode | 39.8 MB | 19.9 MB | 50% |
| + delta-encode phi_sa+run_pos | 19.9 MB | 10.9 MB | 45% more |
| + bit-pack phi_da+run_sa | 10.9 MB | 7.0 MB | 36% more |
| **All combined** | **39.8 MB** | **~7 MB** | **82%** |

For comparison, the SSA at s=8 is 0.3 MB. Even with all optimizations, the SRI would still
be ~50x larger than SSA, because the SRI fundamentally stores O(r) entries while SSA stores
O(n/2^s) entries, and r >> n/256 for non-repetitive data.

## 9. Recommendations

1. **Ship the s=1 dedup and bitvector optimizations** (this worktree). They are safe,
   backward-compatible, and provide 33% disk savings for s=1 with no performance regression.

2. **Implement 32-bit mode as next priority**. It halves index size for n < 4G (covers all
   practical pangenome indexes for years to come). Low-medium complexity.

3. **Consider delta encoding phi_sa/run_pos** for a further ~50% reduction on the sorted
   arrays. This is the biggest remaining win but requires auxiliary index structures.

4. **Do not use s=1** for locate performance. Despite "no subsampling", it walks ~n/r LF
   steps to resolve the toehold, making it 10x slower than s=8 for locate_one. Use s=4 or
   s=8 for the best speed/size tradeoff.

5. **The fundamental size gap** between SRI and SSA comes from storing O(r) entries for the
   phi/toehold structures. For non-repetitive data (r ≈ 0.75n), this is unavoidable. The SRI's
   value proposition is speed, not size - use it when locate throughput justifies the memory cost.
