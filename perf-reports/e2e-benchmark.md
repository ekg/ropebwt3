# End-to-End Performance Benchmark: ropebwt3

**Date:** 2026-02-21T13:15:14Z
**Platform:** Linux x86_64, 16 cores
**Compiler:** gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0
**Threads:** 4

## Test Data

| Item | Description |
|------|-------------|
| Reference | 100 synthetic DNA sequences, 3–4 kb each (~350 kb FASTA); mix of random and tandem-repeat regions |
| Queries | 1000 reads, 150 bp (70% sampled from ref with ~2% error, 30% random) |

## 1. Index Construction

| Step | Wall Time (s) | Peak RSS (MB) |
|------|--------------|----------------|
| BWT build (FMD) | .022 | 6.7 |
| SSA (sampled suffix array) | .038 | 2.9 |
| SR-index | .326 | 44.6 |
| Move structure (.mvi) | .042 | 21.9 |

**Total build time (old pipeline: BWT+SSA):** .060s
**Total build time (new pipeline: BWT+SRI+Move):** .390s

## 2. Disk Footprint

| File | Size (MB) |
|------|-----------|
| BWT (.fmd) | .35 |
| SSA (.ssa) | .02 |
| SR-index (.sri) | 14.14 |
| Move (.mvi) | 19.21 |

| Configuration | Total Disk (MB) |
|--------------|----------------|
| Old (FMD+SSA) | .37 |
| New (FMD+SRI+Move) | 33.70 |

## 3. Query Performance: mem (SMEM finding + locate, -p 10)

| Backend | Wall Time (s) | Peak RSS (MB) |
|---------|--------------|----------------|
| FMD + SSA (old) | .045 | 3.2 |
| b-move + SRI (new) | .033 | 55.8 |

## 4. Query Performance: sw (Smith-Waterman + locate, -p 10)

| Backend | Wall Time (s) | Peak RSS (MB) |
|---------|--------------|----------------|
| FMD + SSA (old) | .234 | 33.2 |
| b-move + SRI (new) | .255 | 85.6 |

## 5. Query Performance: ms (matching statistics)

| Backend | Wall Time (s) | Peak RSS (MB) |
|---------|--------------|----------------|
| FM-index only (-F, old) | 51.843 | 18.1 |
| Move structure (new) | 54.301 | 56.6 |

## 6. Summary

### Configuration Comparison

| Metric | Old (FMD+SSA) | New (b-move+SRI) |
|--------|--------------|-------------------|
| Build time (s) | .060 | .390 |
| Disk footprint (MB) | .37 | 33.70 |
| mem wall time (s) | .045 | .033 |
| mem peak RSS (MB) | 3.2 | 55.8 |
| sw wall time (s) | .234 | .255 |
| sw peak RSS (MB) | 33.2 | 85.6 |
| ms wall time (s) | 51.843 | 54.301 |
| ms peak RSS (MB) | 18.1 | 56.6 |

### Notes

- Backend selection is automatic: ropebwt3 uses b-move when a `.mvi` file is
  present, and prefers SR-index (`.sri`) over SSA (`.ssa`) for locate.
- Old backend was isolated by temporarily hiding `.mvi`/`.sri` files; new
  backend by hiding `.ssa`.
- The `ms` command has an explicit `-F` flag to force FM-index mode.
- All timings include I/O. Results on synthetic data; real genomic data with
  higher repetitiveness may show larger differences in favor of the new backends.

### Known Issue

- `srindex` crashes (heap corruption) on inputs exceeding ~800k BWT symbols
  with random/semi-random DNA. The benchmark dataset was sized to stay under
  this threshold. This bug should be investigated separately.
