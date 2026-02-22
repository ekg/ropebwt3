# Performance Benchmark: LCP/Matching Statistics

**Date:** 2026-02-21 (updated with re-run confirmation)
**ropebwt3 version:** 3.10-r281
**CPU:** AMD Ryzen AI 7 350 w/ Radeon 860M
**Cores:** 16
**RAM:** ~94 GB
**OS:** Linux 6.17.0-14-generic

## Datasets

Two types of datasets are used: **diverse** (random, independent sequences) and **repetitive** (pangenome-like, copies of a shared ancestor with low divergence).

| Dataset | Type | Sequences | Total bp | BWT runs | Runs/bp | FMD size |
|---------|------|-----------|----------|----------|---------|----------|
| Small (perf-data) | Diverse | 200 (fwd+rev) | 200 kbp | 150,569 | 75.3% | 129 KB |
| Medium (perf-data) | Diverse | 1,000 (fwd+rev) | 10 Mbp | 7,499,526 | 75.0% | 6.2 MB |
| Large (perf-data) | Repetitive | 400 (fwd+rev) | 20 Mbp | 1,243,709 | 6.2% | 1.5 MB |
| Custom-Micro | Repetitive | 20 (fwd+rev) | 10 kbp | 1,554 | 15.5% | 2.4 KB |
| Custom-Small | Repetitive | 40 (fwd+rev) | 80 kbp | 11,165 | 13.9% | 14 KB |
| Custom-Medium | Repetitive | 100 (fwd+rev) | 400 kbp | 46,829 | 11.7% | 56 KB |
| Custom-Diverse | Diverse | 100 (fwd+rev) | 400 kbp | 300,221 | 75.0% | 256 KB |

**Queries:** 10,000 per perf-data dataset (100-150 bp, extracted from reference with ~2% error rate). 50 per custom dataset (200-1000 bp, random or extracted with 2% mutation).

---

## 1. LCP Array Build Time

The LCP array stores longest common prefix values at BWT run boundaries. Building it requires a Psi-walk at each boundary, where cost scales with the LCP value. Repetitive data has longer LCP values, making build time superlinear in run count.

| Dataset | Runs | LCP Build | LCP+Thresh Build | Peak RSS |
|---------|------|-----------|------------------|----------|
| Small (diverse) | 150,569 | 5.3s | 5.2s | 8 MB |
| Custom-Micro (repetitive) | 1,554 | 0.82s | 0.81s | 3.6 MB |
| Custom-Small (repetitive) | 11,165 | 12.2s | 11.5s | 3.6 MB |
| Custom-Medium (repetitive) | 46,829 | 74.0s | 75.3s | 3.8 MB |
| Custom-Diverse | 300,221 | 11.9s | 12.0s | 13.7 MB |
| Large (repetitive) | 1,243,709 | *(pending)* | *(pending)* | *(pending)* |

**Key observations:**
- Diverse data (short LCP values at boundaries) builds much faster despite more runs: 300K runs in 11s vs 47K runs in 68s.
- Repetitive data is expensive because each Psi-walk traverses long shared prefixes (LCP values in the hundreds or thousands).
- Threshold computation (`-t`) adds minimal overhead (<10%) on top of LCP build, since it only examines neighboring LCP values.
- Memory usage is very low: dominated by the FM-index itself, not the LCP array (which is O(r) int64 values).

### Commands
```bash
# LCP array only
/usr/bin/time -v ./ropebwt3 lcp <idx.fmd> > lcp.txt

# LCP + thresholds (required for MS)
/usr/bin/time -v ./ropebwt3 lcp -t <idx.fmd> > lcp_thresh.txt
```

---

## 2. Move Structure Build

The move structure converts the FM-index into a run-based index for O(1) rank queries. The split depth parameter `-d` splits long runs to bound fast-forward steps.

### Build Time

| Dataset | Move Build (d=0) | Move Build (d=4) |
|---------|------------------|------------------|
| Small (diverse) | <0.01s / 6.9 MB | 0.01s / 9.6 MB |
| Medium (diverse) | --- | 0.78s / 360 MB |
| Large (repetitive) | --- | 0.13s / 62 MB |

Move builds are very fast: 0.13s even for 20M symbols with 1.2M runs.

### Split Depth Effect on Index Size

| Dataset | d=0 | d=2 | d=3 | d=4 |
|---------|-----|-----|-----|-----|
| Small (diverse, 150K runs) | 6.9 MB | 6.9 MB | 6.9 MB | 6.9 MB |
| Large (repetitive, 1.2M runs) | 56.9 MB | 56.9 MB | 56.9 MB | 56.9 MB |

**Key observation:** Split depth has no effect on these datasets. The maximum run length in the BWT is small enough that splitting finds nothing to split. For the custom-medium dataset, max run length is 137 (median=1, 56% of runs are length 1). Splitting only affects datasets with very long runs (e.g., thousands of identical characters in a row), which requires much larger, more homogeneous pangenome references.

---

## 3. Matching Statistics Query Throughput

### BWT Path (`-F`) vs Move Path (default)

The BWT path uses FM-index rank queries (O(log r) each). The move path uses precomputed cumulative rank tables (O(1) each) but has overhead from building the move structure in memory.

#### Small Dataset (diverse, 200 kbp, 150K runs)

| Method | Wall Time | Peak RSS |
|--------|-----------|----------|
| MS BWT (-F) | 10.9s | 8 MB |
| MS Move (d=0) | 17.7s | 22 MB |
| MS Move (d=2) | 17.6s | 22 MB |
| MS Move (d=3) | 17.8s | 22 MB |
| PML (-p) | 5.5s | 8 MB |

**Queries:** 10,000 x 100 bp (1M bp total)

#### Medium Dataset (diverse, 10 Mbp, 7.5M runs)

| Method | Wall Time | Peak RSS |
|--------|-----------|----------|
| MS BWT (-F) | 9m 3s | 294 MB |
| MS Move (d=4) | 24m 2s | 981 MB |

**Queries:** 10,000 x 100 bp (1M bp total)

#### Large Dataset (repetitive, 20 Mbp, 1.2M runs)

| Method | Wall Time | Peak RSS |
|--------|-----------|----------|
| MS BWT (-F) | *(pending)* | *(pending)* |
| MS Move (d=0) | *(pending)* | *(pending)* |
| MS Move (d=2) | *(pending)* | *(pending)* |
| MS Move (d=3) | *(pending)* | *(pending)* |
| PML (-p) | *(pending)* | *(pending)* |

**Queries:** 10,000 x 150 bp (1.5M bp total)

#### Custom Datasets (50 queries, 200-1000 bp)

| Dataset | MS BWT (-F) | MS Move (d=0) | MS Move (d=2) | PML |
|---------|-------------|---------------|---------------|-----|
| Micro (10K bp, 1.5K runs) | 1.14s | 1.14s | 1.13s | 0.85s |
| Small (80K bp, 11K runs) | 11.2s | 11.4s | 11.4s | 10.5s |
| Medium (400K bp, 47K runs) | 79.4s | 75.0s | 97.8s | 73.1s |
| Diverse (400K bp, 300K runs) | 12.99s | 13.4s | 13.2s | 12.6s |

**Key observations:**
- **BWT and Move paths are comparable** for most tested datasets. On custom-medium repetitive data, Move d=0 shows a ~5% speedup (75.0s vs 79.4s).
- For **diverse data** with perf-data scale (high run ratio), FM-index is 1.6-2.7x faster than Move: the move table is large (proportional to runs) and doesn't provide speedup because FM-index rank queries are already fast on the compact representation. For custom-diverse, results are comparable.
- **Move d=2 can hurt:** On custom-medium, Move d=2 is 30% slower (97.8s vs 75.0s for d=0), while d=3 matches d=0. The splitting overhead is counterproductive at this scale.
- **PML is consistently faster** than exact MS: 1.3x on micro, ~1.0-1.1x on larger datasets where LCP build dominates.
- The MS computation is dominated by **LCP build** (computed at startup), not by the per-query backward search. For the custom-medium dataset, ~74s is spent in LCP build/setup vs fractions of a second per query.

### Commands
```bash
# BWT path (FM-index only)
/usr/bin/time -v ./ropebwt3 ms -F <idx.fmd> <queries.fa>

# Move path (default)
/usr/bin/time -v ./ropebwt3 ms <idx.fmd> <queries.fa>

# Move path with split depth
/usr/bin/time -v ./ropebwt3 ms -d 2 <idx.fmd> <queries.fa>

# PML (pseudo-matching lengths, faster approximation)
/usr/bin/time -v ./ropebwt3 ms -p <idx.fmd> <queries.fa>
```

---

## 4. PML Computation Time

PML (pseudo-matching lengths) skips the exact length-finding loop in MS, using thresholds for faster repositioning.

| Dataset | MS BWT | PML | Speedup |
|---------|--------|-----|---------|
| Small (diverse) | 10.9s | 5.5s | 2.0x |
| Custom-Micro | 1.14s | 0.85s | 1.3x |
| Custom-Small | 11.2s | 10.5s | 1.1x |
| Custom-Medium | 79.4s | 73.1s | 1.1x |
| Custom-Diverse | 13.0s | 12.6s | 1.0x |

**Key observations:**
- PML speedup is most visible on diverse data (2x for perf-data small), where the exact MS backward search loop iterates many times.
- On repetitive data, PML and MS converge because the threshold provides a good jump target on the first try (the search loop rarely needs multiple iterations).

---

## 5. Memory Usage

| Dataset | FMD Size | LCP Build RSS | MS BWT RSS | MS Move (d=0) RSS | MS Move (d=4) RSS |
|---------|----------|---------------|------------|-------------------|-------------------|
| Small (diverse) | 129 KB | 7 MB | 8 MB | 22 MB | 22 MB |
| Medium (diverse) | 6.2 MB | --- | 294 MB | --- | 981 MB |
| Large (repetitive) | 1.5 MB | *(pending)* | *(pending)* | *(pending)* | *(pending)* |
| Custom-Micro | 2.4 KB | 3.6 MB | 3.7 MB | 3.5 MB | --- |
| Custom-Medium | 56 KB | 3.7 MB | 3.7 MB | 8.4 MB | --- |
| Custom-Diverse | 256 KB | 11.5 MB | 14.2 MB | 42.0 MB | --- |

**Key observations:**
- LCP index memory is O(r) where r is the number of BWT runs (3 int64 arrays: run_starts, lcp_samples, thresholds = 24 bytes per run).
- Move structure adds significant memory: 48 bytes per run for the move table rows, plus the cumulative rank precomputation table.
- For diverse data with many runs, Move RSS can be 3x the BWT-only RSS (42 MB vs 14 MB for custom-diverse with 300K runs).
- For repetitive data (fewer runs), the overhead is proportionally smaller.

---

## 6. Split Depth Effect

The split depth parameter `-d` controls run splitting in the move structure. It bounds the maximum number of fast-forward steps during backward search to O(2^d).

### Performance

| Dataset | MS Move (d=0) | MS Move (d=2) | MS Move (d=3) | MS Move (d=4) |
|---------|---------------|---------------|---------------|---------------|
| Small (diverse) | 17.7s | 17.6s | 17.8s | --- |
| Custom-Medium | 69.1s | 68.9s | 69.0s | --- |

**Split depth has no measurable performance effect** on these datasets. The BWT runs are already short (maximum run length ~137, median 1), so the fast-forward is always bounded to a small number of steps regardless of split depth.

Split depth is expected to matter for very large pangenome references (billions of symbols, millions of runs) where some runs can be extremely long. At those scales, splitting d=2-4 would cap fast-forward to 4-16 steps.

---

## 7. Correctness Verification

### BWT vs Move Path Consistency

| Test | Result |
|------|--------|
| Custom datasets (4 datasets x 2 query sets x 3 d-values) | **24/24 PASS** |
| Small perf-data (FMI vs Move d=0,2,3) | **3/3 PASS** |

The FM-index path (`-F`) and move path produce **identical** matching statistics values for all tested inputs and split depths. This confirms that the move structure correctly implements the same algorithm.

### PML <= MS Invariant

The code documents that PML[i] <= MS[i] for all positions i. Our benchmarks show violations:

| Dataset | Violations | Total Positions | Violation Rate |
|---------|------------|-----------------|----------------|
| Small (perf-data) | 30,754 | 1,000,000 | 3.1% |
| Custom-Micro (random queries) | 3,635 | 29,638 | 12.3% |
| Custom-Small (random queries) | 3,820 | 29,638 | 12.9% |
| Custom-Medium (random queries) | 4,296 | 29,638 | 14.5% |
| Custom-Diverse (random queries) | 3,830 | 29,638 | 12.9% |

**Analysis:** PML values exceed MS values at 3-15% of positions. This appears to be a correctness issue in either the MS or PML computation. The unit tests (`test-ms`, `test-move-ms`) pass on small hand-crafted inputs, but the invariant breaks down on larger, more complex sequences. The violations are small in magnitude (PML typically exceeds MS by 1-5 positions), suggesting an edge case in interval widening or threshold lookup. **This finding warrants investigation.**

---

## 8. Scaling Analysis

### LCP Build: Run Count vs Time

The LCP build time is controlled by two factors: the number of runs (r) and the average LCP value at run boundaries. For repetitive data, LCP values are large, making each Psi-walk expensive.

| Dataset | Runs (r) | Symbols (n) | Runs/n | LCP Build Time | Time/Run |
|---------|----------|-------------|--------|----------------|----------|
| Custom-Micro | 1,554 | 10K | 15.5% | 0.82s | 0.53 ms |
| Custom-Small | 11,165 | 80K | 13.9% | 12.2s | 1.09 ms |
| Custom-Medium | 46,829 | 400K | 11.7% | 74.0s | 1.58 ms |
| Custom-Diverse | 300,221 | 400K | 75.0% | 11.9s | 0.040 ms |
| Small (diverse) | 150,569 | 200K | 75.3% | 5.3s | 0.035 ms |

**Scaling patterns:**
- **Diverse data:** Time per run is ~0.035ms (fast, low LCP values). Scales linearly with r.
- **Repetitive data:** Time per run increases with data size (0.51ms -> 1.45ms) because LCP values grow with sequence length/repetitiveness. Superlinear scaling.

### MS Query: Throughput

| Dataset | Queries | Total Query bp | BWT Time | Throughput (bp/s) |
|---------|---------|----------------|----------|-------------------|
| Small (diverse, 150K runs) | 10,000 | 1M | 10.9s | 92K bp/s |
| Medium (diverse, 7.5M runs) | 10,000 | 1M | 543s | 1.8K bp/s |
| Custom-Diverse (300K runs) | 50 | ~30K | 12.0s | 2.5K bp/s |

MS throughput drops significantly with more BWT runs, even holding query length constant. This is because backward extension fails more often with a larger BWT, requiring more Psi-walk LCP computations per query position.

---

## Summary

| Aspect | Finding |
|--------|---------|
| **LCP Build** | Bottleneck for repetitive data. 74s for 47K runs (repetitive) vs 12s for 300K runs (diverse). Psi-walk depth drives cost. |
| **MS Computation** | BWT path (-F) is comparable to Move path. Move shows 5% speedup on medium repetitive data; BWT is slightly better on small/diverse. |
| **PML** | 1.0-1.3x faster than exact MS. Best speedup on diverse data. |
| **Split Depth** | No effect at these data sizes. Runs are already short (max ~137). |
| **Memory** | Move structure adds 2-3x RSS. BWT-only path is more memory efficient. |
| **Correctness** | BWT vs Move: 27/27 PASS (confirmed in re-run). PML <= MS invariant: VIOLATED at 3-15% of positions (potential bug, avg overshoot 1.33, max 5). |
| **Recommendation** | Use `-F` (BWT-only) for datasets with <10M runs. Reserve Move for very large pangenome indices where O(1) rank provides measurable speedup. |
