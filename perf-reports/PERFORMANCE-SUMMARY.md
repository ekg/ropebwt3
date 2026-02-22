# ropebwt3 Performance Summary

**Date:** 2026-02-21
**Version:** 3.10-r281 (commit 9dbf839)
**Platform:** Linux 6.17.0-14-generic, AMD Ryzen AI 7 350, 16 cores, 94 GB RAM

This report synthesizes four component benchmarks (move structure, SR-index, matching statistics, end-to-end) run on synthetic datasets ranging from 10 kbp to 20 Mbp. All findings are on synthetic data; real pangenome-scale datasets may shift the balance toward the new backends.

---

## Overall Verdict

**The new backends (b-move + SR-index) are not yet faster or smaller for typical workloads at the scales tested.** The FMD + SSA pipeline remains the better default for indexes under ~1 GB. However, the SR-index delivers a clear 2-8x speedup for position locate, and the move/b-move structures have sound theoretical properties (O(r)-space, O(1) rank) that should pay off at pangenome scale where indexes far exceed CPU cache.

**Recommendation: keep FMD + SSA as the default pipeline.** Build the SR-index only when MEM locate throughput is the bottleneck. Build the move structure only for very large pangenome indexes where cache effects favor run-length-proportional data structures.

---

## Per-Component Comparison

### Speed

| Operation | Old (FMD+SSA) | New (b-move+SRI) | Ratio | Winner |
|-----------|---------------|-------------------|-------|--------|
| **Index build** | 0.06 s | 0.39 s | 6.5x slower | Old |
| **MEM find (SMEM)** | 0.045 s | 0.033 s | 1.4x faster | New |
| **MEM locate** (low-p) | baseline | 2-2.3x faster | 2x faster | New (SRI) |
| **MEM locate** (high-p) | baseline | 4.6-7.8x faster | 5-8x faster | New (SRI) |
| **SW alignment** | 0.234 s | 0.255 s | ~1x (comparable) | Tie |
| **Matching statistics** | 51.8 s | 54.3 s | ~1x (comparable) | Tie |
| **MS (medium diverse)** | 9 min | 24 min | 2.7x slower | Old |

*MEM find uses b-move for backward search; MEM locate uses SRI's phi-function. The locate speedup is real and grows with the number of reported positions. MS time is dominated by LCP construction, not the query backend.*

### Memory (Peak RSS)

| Operation | Old (FMD+SSA) | New (b-move+SRI) | Ratio |
|-----------|---------------|-------------------|-------|
| MEM query | 3-12 MB | 50-715 MB | 4-60x more |
| SW query | 33 MB | 86-131 MB | 2.6-4x more |
| MS query | 18-301 MB | 57-1004 MB | 3x more |

### Disk

| Component | Size | Notes |
|-----------|------|-------|
| FMD (.fmd) | 0.35-6.2 MB | Run-length delta encoded; very compact |
| SSA (.ssa) | 0.02-0.3 MB | Minimal |
| SRI (.sri) | 14-60 MB | 80-200x larger than SSA |
| MVI (.mvi) | 7-343 MB | 39-55x larger than FMD |
| **Old pipeline total** | **0.37-6.5 MB** | |
| **New pipeline total** | **21-409 MB** | **57-91x larger** |

---

## SR-Index: Recommended Configuration

The SR-index (SRI) replaces the sampled suffix array (SSA) for converting BWT positions to genomic coordinates.

| s value | MEM locate speedup | Index size | Build RSS | Recommendation |
|---------|-------------------|------------|-----------|----------------|
| 1 | 1.3-4.6x | 32 MB | 87 MB | Avoid: slower than s=4/8 despite no subsampling |
| 4 | 2.3x | 60 MB | 164 MB | Fast but large; use if memory is unconstrained |
| **8 (default)** | **2.1-7.8x** | **40 MB** | **114 MB** | **Best balance of speed, size, and memory** |
| 16 | 2.1x | 31 MB | 98 MB | Reasonable if disk-constrained |
| 64 | 1.9x | 24 MB | 90 MB | Smallest SRI; modest speedup |

**Use s=8 (the default).** It provides the best tradeoff: 2-8x faster locate than SSA, with manageable 40 MB index size. The speedup grows with the number of positions requested per query (`-p`), making SRI most valuable for high-coverage MEM searches.

---

## Tradeoffs: When to Use Which Backend

### Use FMD + SSA (old pipeline) when:
- Index fits in L2/L3 cache (FMD < ~1 GB)
- Memory is constrained (old pipeline uses 4-60x less RAM)
- Disk is constrained (old pipeline uses 57-91x less disk)
- Workload is SW alignment (locate is not the bottleneck)
- Workload is matching statistics on diverse data (FMD 1.6-2.7x faster)
- Data is random/low-redundancy (many BWT runs make move table large and cache-hostile)

### Use FMD + SRI (SR-index only, no move) when:
- MEM locate throughput is the bottleneck
- High `-p` values (many positions per hit) amplify the 5-8x speedup
- Memory budget allows ~50 MB RSS overhead for the SRI

### Use b-move + SRI (full new pipeline) when:
- Index exceeds CPU cache capacity (multi-GB FMD)
- Data is highly repetitive (r << n), e.g., large pangenome collections
- O(r)-space guarantee is needed for memory predictability
- Theoretical O(1) rank per LF-step matters (very long backward searches)

### Use PML instead of exact MS when:
- Approximate matching lengths are acceptable
- Diverse data: PML is up to 2x faster than exact MS
- Repetitive data: PML converges to exact MS speed (both dominated by LCP build)

---

## Surprising Findings and Regressions

### 1. PML > MS invariant violation (potential bug)
PML (pseudo-matching lengths) should satisfy PML[i] <= MS[i] by definition, but **3-15% of positions show PML > MS** with an average overshoot of 1.33 and max of 5. This was observed consistently across all datasets and query types. The BWT-vs-Move MS correctness checks all pass (27/27), so the issue is specific to the PML computation or the exact MS computation at edge cases. **This warrants investigation.**

### 2. Split depth has no effect at tested scales
The move structure's split depth parameter (`-d`) had zero measurable impact on performance. BWT runs at these scales are already short (max ~137, median 1). The parameter only matters for very large pangenomes with extremely long runs. Using `-d2` actually caused a 30% regression in one case.

### 3. LCP construction is the MS bottleneck, not the query backend
Matching statistics time is dominated by LCP array construction (a full Psi-walk over all run boundaries), not by the per-query backward search. For repetitive data, LCP build is superlinearly expensive: 47K runs took 74s (1.58 ms/run) vs 300K diverse runs in 12s (0.04 ms/run). This means switching backends (BWT vs Move) barely affects total MS wall time.

### 4. Move/b-move is slower than FMD at all tested scales for SMEM
Despite O(1) theoretical rank cost, b-move was 1.1-4.5x slower than FMD backward search for SMEM finding. The constant factor and cache pressure from the 48-byte-per-row move table outweigh the algorithmic advantage when the FMD fits in cache. The crossover is expected at much larger scales.

### 5. Heap buffer overflow in srindex.c (fixed)
A heap corruption bug was found in `srindex.c:sa_worker_func()` where the per-sentinel target buffer was hardcoded to 4096 entries but could overflow with larger datasets (~11,700 targets/sentinel in the SR-index benchmark). **Fixed** with proportional sizing (4x average + slack). The e2e benchmark initially hit this bug on >800k BWT symbols; the fix resolves it.

---

## Quick Reference

```
Pipeline selection flowchart:

  Is MEM locate the bottleneck?
  ├─ Yes → Build SRI (s=8). Consider b-move only if FMD > 1 GB.
  └─ No
       ├─ SW alignment? → Use FMD + SSA (old pipeline)
       ├─ Matching statistics?
       │   ├─ Exact needed? → Use `ms -F` (BWT path)
       │   └─ Approximate OK? → Use `ms -p` (PML, up to 2x faster)
       └─ Default → FMD + SSA (old pipeline)
```

---

## Source Reports

| Report | File |
|--------|------|
| Move structure benchmark | [perf-reports/move-benchmark.md](move-benchmark.md) |
| SR-index benchmark | [perf-reports/srindex-benchmark.md](srindex-benchmark.md) |
| Matching statistics benchmark | [perf-reports/ms-benchmark.md](ms-benchmark.md) |
| End-to-end benchmark | [perf-reports/e2e-benchmark.md](e2e-benchmark.md) |
