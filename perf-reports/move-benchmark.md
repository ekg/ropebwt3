# Performance Benchmark: Move Structure vs FMD Backward Search

**Date:** 2026-02-21
**ropebwt3 version:** 3.10-r281
**CPU:** AMD Ryzen AI 7 350 w/ Radeon 860M
**Cores:** 16
**RAM:** ~94 GB
**OS:** Linux 6.17.0-14-generic

## Datasets

| Dataset | Description | Sequences | Total bp | BWT runs | FMD size | Compressibility |
|---------|-------------|-----------|----------|----------|----------|-----------------|
| Small   | 100 random seqs x 1 kbp | 200 (fwd+rev) | 200 kbp | ~150k | 129 KB | Low (random) |
| Medium  | 500 random seqs x 10 kbp | 1,000 (fwd+rev) | 10 Mbp | 7,499,526 | 6.2 MB | Low (random) |
| Large   | 200 seqs x 50 kbp (1% divergence from shared ancestor) | 400 (fwd+rev) | 20 Mbp | 1,243,709 | 1.5 MB | High (redundant) |

**Queries:** 10,000 per dataset, 100-150 bp each, extracted from reference with ~2% error rate.

---

## 1. Index Build Time

| Dataset | FMD Build | Move Build (d=4) | Move Overhead |
|---------|-----------|-------------------|---------------|
| Small   | 0.01s / 3.7 MB | 0.01s / 9.6 MB | ~0s |
| Medium  | 0.24s / 51 MB | 0.78s / 360 MB | +0.54s (3.3x) |
| Large   | 0.40s / 101 MB | 0.13s / 62 MB | -0.27s (0.3x) |

**Key observations:**
- FMD build (via `ropebwt3 build -d`) includes BWT construction from FASTA using libsais.
- Move index build (`ropebwt3 move -d4`) converts the FMD to a move table with run-splitting depth 4.
- For the redundant large dataset, move build is faster because there are far fewer runs (1.24M vs 7.5M for medium).
- Move build peak RSS is proportional to the number of runs (each row is 48 bytes + run splitting).

### Commands
```bash
# Build FMD
/usr/bin/time -v ./ropebwt3 build -do perf-data/small.fmd perf-data/small.fa

# Build Move index from FMD (split depth 4)
/usr/bin/time -v ./ropebwt3 move -d4 perf-data/small.fmd perf-data/small.fmd.mvi
```

---

## 2. Index Size on Disk

| Dataset | FMD (.fmd) | MVI (.fmd.mvi) | MVI/FMD Ratio | BWT Runs |
|---------|-----------|-----------------|---------------|----------|
| Small   | 129 KB    | 6.9 MB          | 55x           | ~150k    |
| Medium  | 6.2 MB    | 343 MB          | 55x           | 7.5M     |
| Large   | 1.5 MB    | 56.9 MB         | 39x           | 1.24M    |

**Key observations:**
- The MVI format is 39-55x larger than FMD. Each move row is 48 bytes, and with split depth 4, runs are subdivided further.
- FMD is run-length delta encoded — extremely compact for redundant data (large dataset: 1.5 MB for 20 Mbp).
- The move table's size is proportional to the number of (split) runs × 48 bytes, plus overhead.
- For large pangenome datasets where r << n, the MVI size relative to BWT length is reasonable, but absolute size remains much larger than FMD.

---

## 3. SMEM Finding Throughput (`ropebwt3 mem`)

The `mem` command automatically uses b-move for SMEM finding when a `.mvi` file is present alongside the `.fmd`. Otherwise it uses the FMD backend directly.

### Single-threaded (10,000 queries)

| Dataset | FMD (wall) | b-move (wall) | Ratio | FMD RSS | b-move RSS |
|---------|-----------|---------------|-------|---------|------------|
| Small   | 0.21s     | 0.27s         | 1.3x slower | 5.9 MB | 19.6 MB |
| Medium  | 0.29s     | 1.02s         | 3.5x slower | 12 MB  | 715 MB  |
| Large   | 0.44s     | 1.02s         | 2.3x slower | 7.9 MB | 124 MB  |

### Multi-threaded, 4 threads (10,000 queries)

| Dataset | FMD (wall) | b-move (wall) | Ratio | FMD RSS | b-move RSS |
|---------|-----------|---------------|-------|---------|------------|
| Small   | 0.06s     | 0.08s         | 1.3x slower | 5.7 MB | 19.7 MB |
| Medium  | 0.08s     | 0.36s         | 4.5x slower | 11.9 MB | 715 MB |
| Large   | 0.14s     | 0.30s         | 2.1x slower | 7.7 MB | 124 MB |

**Key observations:**
- FMD backward search is consistently faster than b-move for SMEM finding on all dataset sizes.
- The gap is largest on the medium (random) dataset: b-move is 3.5-4.5x slower. This is because random data produces many BWT runs, making the move table large and causing more cache misses.
- For redundant data (large), the gap narrows to ~2x because the move table is smaller and more cache-friendly.
- Both produce identical SMEM results (same hit counts).
- Memory overhead of b-move is significant: the entire MVI must be loaded into memory (or mmapped).

### Commands
```bash
# SMEM with FMD (no .mvi file present)
/usr/bin/time -v ./ropebwt3 mem -t1 perf-data/small.fmd perf-data/small_queries.fa

# SMEM with b-move (.mvi file alongside .fmd)
# First: cp perf-data/small.mvi.bak perf-data/small.fmd.mvi
/usr/bin/time -v ./ropebwt3 mem -t1 perf-data/small.fmd perf-data/small_queries.fa
```

---

## 4. Count Query Throughput (`ropebwt3 mem -l1 -c1`)

Using `mem -l1 -c1` tests backward search with minimal length filter, stressing the core rank/LF-mapping.

### Single-threaded (10,000 queries, 100 bp each)

| Dataset | FMD (wall) | b-move (wall) | Ratio | FMD RSS | b-move RSS |
|---------|-----------|---------------|-------|---------|------------|
| Small   | 0.41s     | 0.45s         | 1.1x slower | 9.1 MB | 23.4 MB |
| Medium  | 0.70s     | 1.83s         | 2.6x slower | 16 MB  | 719 MB  |
| Large   | (not measured separately — see SMEM results for comparable data) | | | | |

**Key observations:**
- With `-l1 -c1`, more matches are reported, increasing output overhead.
- The FMD backend is faster for pure count-style queries, with the gap widening on larger, random data.
- For the small dataset the performance difference is minimal (~10%).

### Commands
```bash
/usr/bin/time -v ./ropebwt3 mem -t1 -l1 -c1 perf-data/small.fmd perf-data/small_queries.fa
```

---

## 5. Matching Statistics (`ropebwt3 ms`)

The `ms` command computes matching statistics using LCP array thresholds. With `-F`, only the FMD is used for backward search during LCP construction and MS computation. Without `-F`, a move table is built on-the-fly and used for the backward search steps.

**Note:** The LCP array construction is the dominant cost. It requires O(n) backward search steps through the BWT regardless of backend, where n is the BWT length. The actual query phase (computing MS for input reads) is fast in comparison.

### Results (10,000 queries, 100 bp each)

| Dataset | FMD `-F` (wall) | Move `-d4` (wall) | Ratio | FMD RSS | Move RSS |
|---------|-----------------|-------------------|-------|---------|----------|
| Small   | 10.88s          | 18.27s            | 1.7x slower | 8.0 MB | 22.3 MB |
| Medium  | 9:03 (543s)     | 24:02 (1442s)     | 2.7x slower | 301 MB | 1004 MB |
| Large   | >37 min (killed)| (skipped)         | — | — | — |

**Key observations:**
- The FMD backend is 1.7-2.7x faster for matching statistics computation.
- The time is dominated by LCP construction (backward search through the entire BWT), not query processing.
- The large redundant dataset was killed after 37+ minutes. Despite fewer runs, the 20M BWT symbols require proportionally more backward search steps. The high redundancy causes very long LCP values, increasing total work.
- Memory: move backend uses ~3x more RSS due to the on-the-fly move table construction plus LCP arrays.

### Commands
```bash
# MS with FMD backend
/usr/bin/time -v ./ropebwt3 ms -F perf-data/medium.fmd perf-data/medium_queries.fa

# MS with move backend (split depth 4)
/usr/bin/time -v ./ropebwt3 ms -d4 perf-data/medium.fmd perf-data/medium_queries.fa
```

---

## 6. Memory Usage Summary (Peak RSS)

| Operation | Dataset | FMD | Move/b-move |
|-----------|---------|-----|-------------|
| Index build | Small | 3.7 MB | 9.6 MB |
| Index build | Medium | 52 MB | 360 MB |
| Index build | Large | 101 MB | 62 MB |
| SMEM (1 thread) | Small | 5.9 MB | 19.6 MB |
| SMEM (1 thread) | Medium | 12 MB | 715 MB |
| SMEM (1 thread) | Large | 7.9 MB | 124 MB |
| MS | Small | 8.0 MB | 22.3 MB |
| MS | Medium | 301 MB | 1004 MB |

**Key observations:**
- FMD is extremely memory-efficient due to run-length delta encoding.
- Move/b-move requires loading or building the full move table (48 bytes × n_split_runs).
- For the redundant large dataset, move memory is moderate (124 MB for SMEM) because there are only 1.24M runs.
- For random data (medium), move memory is high (715 MB) due to 7.5M runs × split factor.

---

## 7. Analysis and Conclusions

### When to prefer FMD (current default)
- **Random or low-redundancy data:** FMD is significantly faster and much more memory-efficient.
- **Memory-constrained environments:** FMD indexes are 40-55x smaller than move indexes.
- **Standard SMEM/count queries:** FMD backward search is 1.1-4.5x faster than b-move across all tested scenarios.

### When move structure may be advantageous
- **Highly redundant data at scale:** The move table scales with BWT runs (r), not BWT length (n). For pangenome-scale datasets (millions of sequences, r << n), the move table may become more competitive.
- **Cache-friendly sequential access:** The move structure is designed for sequential LF-mapping in run-order, which can be more cache-friendly on very large indexes that don't fit in LLC. This advantage was not observed at our test scale but is expected on multi-gigabyte indexes.
- **Theoretical guarantee:** Move provides O(r)-space backward search with O(1) LF-mapping per step (after repositioning), vs. FMD's O(n/block_size) space with O(block_decode_cost) per rank query. The crossover point where move becomes faster is at much larger scales than tested.

### Limitations of this benchmark
1. **Synthetic data only:** Real genomic data (e.g., 152 M. tuberculosis genomes or 579 human assemblies) would stress different aspects, particularly cache behavior on multi-GB indexes.
2. **Small scale:** At 5-20 Mbp, the FMD fits entirely in L2/L3 cache, giving it a major advantage. The move structure's benefits emerge when the index exceeds cache capacity.
3. **LCP construction dominates MS benchmarks:** The matching statistics time is dominated by LCP construction, not query processing. A pre-built LCP index would change these results dramatically.
4. **No pre-built MVI for MS:** The `ms` command builds the move table on-the-fly from FMD. Using a pre-built `.mvi` (if supported in the future) would reduce startup cost.

### Recommendations for production use
- For indexes under ~1 GB FMD size: **use FMD only** (do not build .mvi)
- For very large pangenome indexes (>10 GB FMD): benchmark with the actual data, as move/b-move may become competitive due to cache effects
- The move structure is still valuable for the `ms` command's theoretical guarantees (O(r)-space matching statistics via MONI algorithm), even if the constant factor is currently higher

---

## 8. Reproducibility

### Generate test data
```bash
python3 perf-data/gen_data.py
```

### Run full benchmark
```bash
bash perf-data/run_benchmark.sh 2>&1 | tee perf-reports/benchmark.log
```

### Run individual tests
```bash
# Build FMD index
./ropebwt3 build -do perf-data/medium.fmd perf-data/medium.fa

# Build move index (split depth 4)
./ropebwt3 move -d4 perf-data/medium.fmd perf-data/medium.fmd.mvi

# SMEM with FMD (ensure no .mvi file alongside .fmd)
./ropebwt3 mem -t1 perf-data/medium.fmd perf-data/medium_queries.fa

# SMEM with b-move (ensure .mvi file alongside .fmd)
./ropebwt3 mem -t1 perf-data/medium.fmd perf-data/medium_queries.fa

# Matching statistics (FMD)
./ropebwt3 ms -F perf-data/medium.fmd perf-data/medium_queries.fa

# Matching statistics (move, split depth 4)
./ropebwt3 ms -d4 perf-data/medium.fmd perf-data/medium_queries.fa

# BWT statistics
./ropebwt3 stat perf-data/medium.fmd
```

### Dataset characteristics
```
$ ./ropebwt3 stat perf-data/small.fmd    # 200 seqs, 200k symbols, ~150k runs
$ ./ropebwt3 stat perf-data/medium.fmd   # 1000 seqs, 10M symbols, 7.5M runs
$ ./ropebwt3 stat perf-data/large.fmd    # 400 seqs, 20M symbols, 1.24M runs
```
