# SR-Index vs Sampled Suffix Array: Performance Benchmark

**Date:** 2026-02-21
**ropebwt3 version:** 3.10-r281 (commit 9dbf839, with srindex buffer fix)
**Platform:** Linux 6.17.0-14-generic, 4 threads

## Test Dataset

| Property          | Value                                                    |
|-------------------|----------------------------------------------------------|
| Reference         | 5 chromosomes x 20 haplotypes (0.5% divergence), 5 Mbp  |
| BWT length (n)    | 10,000,200 symbols                                      |
| BWT runs (r)      | 700,020                                                  |
| Sequences (m)     | 200 (100 forward + 100 reverse complement)               |
| r/n ratio         | 0.070 (moderately repetitive)                            |
| Queries           | 5,000 sequences, 50-500 bp, 70% exact / 30% mutated     |
| FMD index size    | 1.3 MB                                                   |

## Bug Fix During Benchmarking

A heap buffer overflow was discovered and fixed in `srindex.c` during this benchmarking work.
The per-sentinel target buffer in `sa_worker_func` was capped at 4096 entries, but with
larger datasets each sentinel walk can hit many more run-boundary targets. The fix uses a
proportional estimate (4x average per sentinel + slack, capped at total targets).

## 1. Index Build Time

All builds use 4 threads. Wall-clock times are median of 3 runs.

| Index     | s parameter | Wall time (s) | CPU time (s) | Peak RSS (MB) |
|-----------|-------------|---------------|--------------|----------------|
| **SSA**   | 8           | 0.62          | 2.46         | 5.6            |
| **SRI**   | 1           | 2.10          | 7.10         | 87.0           |
| **SRI**   | 4           | 2.49          | 7.53         | 164.1          |
| **SRI**   | 8 (default) | 2.42          | 7.77         | 114.0          |
| **SRI**   | 16          | 2.25          | 7.43         | 97.5           |
| **SRI**   | 64          | 2.16          | 7.26         | 89.9           |

**Analysis:** SSA build is ~3.5x faster than any SRI variant. SRI build time is dominated by
the backward walk from each sentinel to compute SA values at run boundaries, which is
fundamentally O(n) work regardless of `s`. The `s` parameter primarily affects memory during
build (collecting subsampled positions) and has minimal impact on build time. Peak RSS during
build increases significantly for small `s` values (more subsampled entries to store).

## 2. Index Size on Disk

| Index     | s parameter | Size (bytes)  | Size (human) | Relative to FMD |
|-----------|-------------|---------------|--------------|-----------------|
| **FMD**   | -           | 1,328,952     | 1.3 MB       | 1.0x            |
| **SSA**   | 8           | 314,132       | 307 KB       | 0.24x           |
| **SRI**   | 1           | 33,604,216    | 32 MB        | 25.3x           |
| **SRI**   | 4           | 62,407,096    | 60 MB        | 47.0x           |
| **SRI**   | 8 (default) | 42,407,096    | 40 MB        | 31.9x           |
| **SRI**   | 16          | 32,407,096    | 31 MB        | 24.4x           |
| **SRI**   | 64          | 24,906,296    | 24 MB        | 18.7x           |

**SRI size breakdown:** The SRI stores:
- Phi arrays: 2 x r x 8 bytes = ~11.2 MB (fixed, for phi_sa and phi_da)
- Run boundary SA: 2 x n_samples x 8 bytes (n_samples ~ 2r = ~1.4M entries)
- Subsampled positions: 2 x n_sub x 8 bytes (n_sub ~ n/s for s > 1)
- Metadata: cumulative lengths, sentinel positions

The `s=4` index is largest because its subsampled set (n/4 ~ 2.5M entries) dominates.
As `s` increases, fewer subsamples are stored, reducing size. At `s=1`, no separate subsampled
array is needed (run boundary samples are reused), making it smaller than `s=4`.

**Key finding:** SSA is 80-200x smaller than SRI for this dataset. The SRI's O(r) base
storage is substantial when r is large relative to n/2^s.

## 3. MEM Locate Performance (mem -p)

### 3.1 Low position count (mem -p5)

Wall-clock time, median of 3 runs, 4 threads, 5000 queries.

| Index     | s   | Wall time (s) | CPU time (s) | Peak RSS (MB) | Speedup vs SSA |
|-----------|-----|---------------|--------------|----------------|----------------|
| **SSA**   | 8   | 0.23          | 0.92         | 8.1            | 1.0x           |
| **SRI**   | 1   | 0.18          | 0.68         | 40.7           | 1.3x           |
| **SRI**   | 4   | 0.10          | 0.34         | 69.0           | 2.3x           |
| **SRI**   | 8   | 0.11          | 0.37         | 49.4           | 2.1x           |
| **SRI**   | 16  | 0.11          | 0.38         | 39.6           | 2.1x           |
| **SRI**   | 64  | 0.12          | 0.44         | 32.0           | 1.9x           |

### 3.2 High position count (mem -p50)

| Index     | s   | Wall time (s) | CPU time (s) | Peak RSS (MB) | Speedup vs SSA |
|-----------|-----|---------------|--------------|----------------|----------------|
| **SSA**   | 8   | 0.78          | 3.09         | 8.5            | 1.0x           |
| **SRI**   | 1   | 0.17          | 0.63         | 41.2           | 4.6x           |
| **SRI**   | 8   | 0.10          | 0.34         | 49.9           | 7.8x           |

**Analysis:** SRI is consistently faster than SSA for MEM locate operations:
- At `-p5`, SRI achieves 1.3-2.3x speedup over SSA.
- At `-p50`, SRI's advantage grows to **4.6-7.8x** because SSA must perform O(n/2^s) LF-mapping
  steps per locate while SRI uses the phi function with O(log r) binary search.
- SRI `s=4` is fastest at `-p5` (sweet spot for this dataset), while `s=8` is fastest at `-p50`.
- The `s=1` variant (no subsampling) is slower than `s=4`/`s=8` because it relies solely on
  phi walks without the shortcut of nearby subsampled positions.
- Trade-off: SRI uses 4-8x more memory at query time due to loading the full index.

## 4. SW Locate Performance (sw -p5)

Wall-clock time, median of 3 runs, 4 threads, 5000 queries.

| Index     | s   | Wall time (s) | CPU time (s) | Peak RSS (MB) | Speedup vs SSA |
|-----------|-----|---------------|--------------|----------------|----------------|
| **SSA**   | 8   | 2.15          | 8.55         | 69.9           | 1.0x           |
| **SRI**   | 1   | 2.16          | 8.53         | 102.6          | 1.0x           |
| **SRI**   | 4   | 2.04          | 8.01         | 130.6          | 1.1x           |
| **SRI**   | 8   | 2.03          | 8.01         | 111.0          | 1.1x           |
| **SRI**   | 16  | 2.12          | 8.35         | 101.2          | 1.0x           |
| **SRI**   | 64  | 2.09          | 8.22         | 93.9           | 1.0x           |

**Analysis:** For SW alignment, locate is a tiny fraction of total work (dominated by the
Smith-Waterman DAWG traversal). SRI provides negligible speedup (<5%) because position lookup
is not the bottleneck.

## 5. Memory Usage During Queries

| Index     | s   | MEM -p5 RSS (MB) | MEM -p50 RSS (MB) | SW -p5 RSS (MB) |
|-----------|-----|-------------------|--------------------|--------------------|
| **SSA**   | 8   | 8.1               | 8.5                | 69.9               |
| **SRI**   | 1   | 40.7              | 41.2               | 102.6              |
| **SRI**   | 4   | 69.0              | -                  | 130.6              |
| **SRI**   | 8   | 49.4              | 49.9               | 111.0              |
| **SRI**   | 16  | 39.6              | -                  | 101.2              |
| **SRI**   | 64  | 32.0              | -                  | 93.9               |

**Analysis:** SSA is significantly more memory-efficient during queries. The SRI must load its
entire phi/run/subsample arrays into RAM, while SSA needs only ~300 KB for its compact sample
table. For MEM queries, SRI uses 4-8.5x more memory; for SW, the overhead is 1.3-1.9x
(smaller relative increase because SW itself uses more working memory).

## 6. Correctness Verification

SSA and SRI produce **identical** results:
- Occurrence counts match exactly for all 5,000 queries
- Position coordinates are equivalent (same genomic locations)
- Minor differences in output: when multiple haplotypes match, SSA and SRI may enumerate
  positions in different order and may report the same location on opposite strands
  (e.g., `chr1_hap2:+:789` vs `chr1_hap5:-:49116` with seqlen=50000 and match_len=95:
  50000 - 789 - 95 = 49116)

## 7. Space/Time Tradeoff Summary

```
  Query Speed (MEM -p5)           Index Size
  ========================        ========================
  SRI s=4  ████████████ 2.3x      SSA s=8  █ 307 KB
  SRI s=8  ██████████   2.1x      SRI s=64 █████████████ 24 MB
  SRI s=16 ██████████   2.1x      SRI s=16 ████████████████ 31 MB
  SRI s=64 █████████    1.9x      SRI s=1  █████████████████ 32 MB
  SRI s=1  ██████       1.3x      SRI s=8  █████████████████████ 40 MB
  SSA s=8  █████        1.0x      SRI s=4  ██████████████████████████████ 60 MB
```

## 8. Recommendations

1. **For speed-critical MEM locate:** Use SRI with `s=4` or `s=8`. Delivers 2-8x speedup
   over SSA for position lookup, with the advantage growing as more positions are requested.

2. **For space-constrained environments:** Use SSA. It is 80-200x smaller than SRI and
   adequate for workloads where position lookup is not the bottleneck (e.g., SW alignment).

3. **Avoid `s=1`:** Despite being the theoretically optimal r-index (no subsampling), it is
   slower than `s=4`/`s=8` for locate because it must take more phi steps without the
   shortcut of nearby subsampled SA values.

4. **SW alignment users:** Stick with SSA unless memory is not a concern. The ~5% speedup
   from SRI does not justify the 60-200x increase in index size.

5. **The sweet spot** for this dataset (r/n ~ 0.07) is SRI `s=8` (the default): good balance
   of 2.1x MEM speedup, 40 MB index size, and 50 MB query RSS. For more repetitive datasets
   (lower r/n), the SRI advantage should increase further as the phi-based lookup becomes
   more efficient.

## 9. Bug Found

During benchmarking, a **heap buffer overflow** was found in `srindex.c:sa_worker_func()`.
The per-sentinel target buffer was hardcoded to `min(n_targets, 4096)` entries, but with
datasets having many run boundaries, a single sentinel's backward walk can hit far more than
4096 target positions (in this dataset, ~11,700 per sentinel on average). This caused heap
corruption and crashes on any non-trivial input.

**Fix:** Changed buffer sizing to `min(n_targets, n_targets / n_sent * 4 + 1024)`, using a
proportional estimate with 4x margin and slack. All existing tests continue to pass.
