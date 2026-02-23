# Optimization Integration Report

**Date:** 2026-02-22
**Base commit:** 4020a55 (master)
**Integration commit:** efd85e5
**Platform:** Linux 6.17.0-14-generic, AMD Ryzen AI 7 350, 16 cores, 94 GB RAM

## Summary

Two study tasks investigated optimization opportunities:

1. **b-move study** (`perf-reports/bmove-analysis.md`): Analyzed the 48-byte/row move table structure, cumrank table memory overhead, and binary search cache miss bottleneck. **No code changes** -- pure analysis identifying future optimization targets (row packing, sampled cumrank, separate p[] array).

2. **SR-index study** (`perf-reports/srindex-analysis.md`): Analyzed the SRI file format inefficiencies and implemented two optimizations:
   - Eliminate s=1 sub array duplication (sub_pos/sub_sa alias run_pos/run_sa)
   - Add bitvector for O(1) sub_pos membership test in `locate_one`

Only the SR-index code changes (commit f08c931) were cherry-picked to master, as the b-move study produced analysis only.

## Changes Applied

### SR-index: s=1 sub array dedup (srindex.c, srindex.h)

**Problem:** For subsampling parameter s=1, `build_subsampled()` allocated new arrays and memcpy'd run_pos/run_sa into sub_pos/sub_sa. This wasted 33% of disk space and memory since the arrays are identical.

**Fix:** For s=1, sub_pos/sub_sa now alias run_pos/run_sa (shared pointers). The serialization format is bumped to v2 (magic `SRI\2`) which writes n_sub=0 for s=1 to indicate the alias. Backward compatibility with v1 files is maintained.

### SR-index: bitvector for locate_one (srindex.c, srindex.h)

**Problem:** `locate_one` performed O(log n_sub) binary search at every LF step to check if the current BWT position is in sub_pos[]. For s=8 with n_sub=2.5M, that's ~21 comparisons per step.

**Fix:** Build a bitvector of n bits (n/8 bytes) marking all sub_pos positions. In `locate_one`, check the bitvector first (O(1) single-word AND), and only do the binary search on a hit.

## Conflict Resolution

No conflicts. The b-move study produced no code changes, and the SR-index changes touch only `srindex.c` and `srindex.h`, which are independent of move.c/move.h.

## Test Results

All 5 test suites pass on the integrated build:

| Test | Result |
|------|--------|
| test-move | PASS (all move/b-move/SMEM tests) |
| test-srindex | PASS (all 28 unit tests) |
| test-ms | PASS (matching statistics) |
| test-move-ms | PASS (move-based matching statistics) |
| test-lcp-threshold | PASS (LCP threshold computation) |

Backward compatibility verified: the new binary correctly reads v1 SRI files (magic `SRI\1`).

## Benchmark Results

### SR-index Disk Size (before vs after)

| Dataset | s | Before (v1) | After (v2) | Change |
|---------|---|-------------|------------|--------|
| Large (20 Mbp, r=1.24M) | 1 | 59,704,488 (56.9 MB) | 39,805,144 (38.0 MB) | **-33.3%** |
| Large (20 Mbp, r=1.24M) | 8 | 79,811,544 (76.1 MB) | 79,811,544 (76.1 MB) | 0% (expected) |
| Small (200 kbp, r=150K) | 1 | 7,230,568 (6.9 MB) | 4,821,464 (4.6 MB) | **-33.3%** |
| Small (200 kbp, r=150K) | 8 | 5,224,664 (5.0 MB) | 5,224,664 (5.0 MB) | 0% (expected) |

The s=1 dedup saves exactly 33% for s=1 indexes (eliminates the redundant sub_pos/sub_sa copy). No change for s>1, as expected.

### SR-index Build Time (large dataset, single-threaded, 3-run median)

| Configuration | Before | After | Change |
|---------------|--------|-------|--------|
| s=1, single-threaded | 10.27s | 9.87s | ~-4% (within noise) |
| s=8, single-threaded | 10.34s | 10.22s | ~-1% (within noise) |

Build time is essentially unchanged. The optimization affects serialization (smaller file) and query time (bitvector), not the SA computation that dominates build.

### SR-index Build Memory

| Configuration | Before RSS | After RSS | Change |
|---------------|-----------|----------|--------|
| s=1, single-threaded | 109,672 KB | 109,520 KB | -0.1% (negligible) |
| s=8, single-threaded | 168,312 KB | 170,652 KB | +1.4% (bitvector overhead, within noise) |

Memory usage is essentially unchanged. The bitvector adds n/8 bytes (2.5 MB for n=20M), but this is small compared to the overall working set.

### locate_one Performance (bitvector optimization)

The bitvector optimization replaces O(log n_sub) binary search per LF step with O(1) bitvector check. This is a constant-factor improvement that is difficult to measure with the current test infrastructure (test-srindex completes in <1ms on small inputs). The theoretical speedup per `locate_one` call:

- **Before:** Each LF step costs O(log n_sub) comparisons (~21 for n_sub=2.5M)
- **After:** Each LF step costs O(1) bitvector check (1 memory access + bitmask)
- **Expected speedup:** ~10-20x per `locate_one` step for large n_sub values

This improvement would be visible on real workloads (e.g., whole-genome locate with millions of positions) but is not measurable with the current benchmark datasets.

### SMEM/FMD/Move Performance (no change expected)

SMEM search output is bit-identical before and after the integration. The SR-index changes are confined to `srindex.c/h` and do not affect the FMD, move, or SMEM code paths.

## Remaining Issues and Future Work

### b-move optimizations (NOT implemented -- analysis only)

The b-move study identified three high-impact optimizations that were NOT implemented:

1. **Pack move rows from 48 to 16 bytes** (HIGH IMPACT, MEDIUM EFFORT): Separate `p[]` into a contiguous array, pack remaining fields. Would reduce .mvi files by 3x and improve cache utilization. Requires new .mvi format version.

2. **Sample cumrank table** (HIGH IMPACT, LOW EFFORT): Store cumrank every 64 rows instead of every row. Would reduce cumrank memory from 343 MB to 5 MB for the medium dataset. Trade-off: rank queries scan ~64 rows sequentially instead of pure lookup.

3. **Separate p[] array for binary search** (MEDIUM IMPACT, LOW EFFORT): Binary search over stride-8 p[] instead of stride-48 rows[]. Would touch 6x fewer cache lines per search.

These are documented in `perf-reports/bmove-analysis.md` for future implementation.

### SR-index optimizations (partially implemented)

The SR-index study identified additional optimizations NOT implemented:

1. **32-bit integer mode** (HIGH IMPACT: ~50% size reduction): For n < 2^32, all arrays could use uint32_t instead of int64_t. Would halve index size.

2. **Delta encoding phi_sa/run_pos** (HIGH IMPACT on sorted arrays): ~19x compression on sorted arrays.

3. **Bit-packed phi_da/run_sa** (MEDIUM IMPACT): Use ceil(log2(n)) bits per entry.

These are documented in `perf-reports/srindex-analysis.md`.

## Conclusion

The integration applies a **safe, well-tested optimization** to the SR-index that achieves a **real 33% disk size reduction for s=1 indexes** with no performance regression. The bitvector optimization provides a theoretical speedup for `locate_one` that will matter on real workloads.

The b-move optimizations require more invasive changes (new file format, restructured data layout) and were correctly kept as analysis-only for future implementation.

**Decision: Push to master.** The SR-index optimizations are:
- Correct (all 28 unit tests pass)
- Backward-compatible (reads v1 SRI files)
- Measurably beneficial (33% smaller s=1 files)
- Non-regressive (no change to SMEM/FMD/move performance)
