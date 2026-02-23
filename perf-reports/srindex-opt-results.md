# SR-Index V3 Compression Results

**Date:** 2026-02-23
**Commit:** worktree `.claude/worktrees/srindex-opt` based on 1d9a3b5

## Summary

Implemented V3 compressed SR-index format with three optimizations:

1. **Delta-encoded sorted arrays** (phi_sa, run_pos, sub_pos): sampled absolute uint32 every 64 entries + 16-bit deltas between samples
2. **Bit-packed unsorted arrays** (phi_da, run_sa, sub_sa): ceil(log2(n)) bits per entry
3. **Compact header** with metadata (bit_width, delta_bits)

Backward compatibility: V3 reader still loads V1/V2 files.

## Size Comparison

### Large dataset (repetitive, n=20,000,400, r=1,243,709, r/n=0.062)

| s | V2 (MB) | V3 (MB) | Reduction |
|---|---------|---------|-----------|
| 1 | 37.9 | 12.4 | **67.3%** |
| 4 | 114.2 | 37.1 | **67.5%** |
| 8 | 76.1 | 24.8 | **67.4%** |
| 16 | 57.0 | 18.6 | **67.4%** |
| 64 | 42.7 | 14.0 | **67.2%** |

### Small dataset (non-repetitive, n=200,200, r=150,569, r/n=0.752)

| s | V2 (KB) | V3 (KB) | Reduction |
|---|---------|---------|-----------|
| 1 | 4,708 | 1,290 | **72.6%** |
| 4 | 5,493 | 1,501 | **72.7%** |
| 8 | 5,102 | 1,396 | **72.6%** |
| 16 | 4,905 | 1,343 | **72.6%** |
| 64 | 4,758 | 1,303 | **72.6%** |

## Load Time Comparison (large dataset)

| s | V2 (ms) | V3 (ms) | Overhead |
|---|---------|---------|----------|
| 1 | 8.3 | 31.4 | 3.8x |
| 8 | 17.9 | 65.7 | 3.7x |
| 64 | 8.8 | 34.7 | 3.9x |

Load time overhead is a one-time cost per index load, well under 100ms for all configurations.

## V3 Format Details

```
Header (52 bytes):
  magic "SRI\3"      4 bytes
  s                   4 bytes (int32)
  m                   8 bytes (int64)
  n                   8 bytes (int64)
  n_runs              8 bytes (int64)
  n_samples           8 bytes (int64)
  n_sub               8 bytes (int64, 0 for s=1 alias)
  bit_width           1 byte  (ceil(log2(n)))
  delta_bits          1 byte  (16 or 32)
  reserved            2 bytes

Arrays:
  phi_sa     delta-encoded (sampled uint32 every 64 + uint16 deltas)
  phi_da     bit-packed (bit_width+1 bits, handles -1 sentinel)
  run_pos    delta-encoded
  run_sa     bit-packed (bit_width bits)
  sub_pos    delta-encoded (if not alias)
  sub_sa     bit-packed (if not alias)
  cum_len    raw int64 (small, m+1 entries)
  tosid      raw int64 (small, m entries)
```

### Delta encoding scheme
- Absolute uint32 sample every K=64 entries
- Between samples: uint16 deltas (or uint32 if any delta > 65535)
- Binary search: O(log(count/K)) on samples + O(K) scan within block

### Bit-packing scheme
- Each value stored in ceil(log2(n)) bits
- Packed into byte array, accessed via bit shifting
- phi_da uses bit_width+1 bits to accommodate -1 sentinel value

## Compression Breakdown

For the large dataset (s=8), V2 stores 79.8 MB of raw int64 arrays:

| Component | V2 bytes | V3 bytes | Savings |
|-----------|----------|----------|---------|
| phi_sa (sorted, 1.24M) | 9,949,672 | ~2,560,000 | 74% (delta) |
| phi_da (unsorted, 1.24M) | 9,949,672 | ~4,040,000 | 59% (bit-pack 26b) |
| run_pos (sorted, 1.24M) | 9,949,672 | ~2,560,000 | 74% (delta) |
| run_sa (unsorted, 1.24M) | 9,949,672 | ~3,880,000 | 61% (bit-pack 25b) |
| sub_pos (sorted, 2.50M) | 20,003,200 | ~5,140,000 | 74% (delta) |
| sub_sa (unsorted, 2.50M) | 20,003,200 | ~7,810,000 | 61% (bit-pack 25b) |
| cum_len + tosid | 6,456 | 6,456 | 0% (raw) |

## Correctness

- All 28 test configurations pass (7 strings × 4 s-values)
- Serialization round-trip verified for every test (dump + restore + verify all arrays + locate_all)
- V1/V2 backward compatibility preserved

## vs. Projected (from srindex-analysis.md)

The analysis projected 82% reduction with all optimizations including 32-bit integer mode. Our V3 achieves 67-73% via delta-encoding and bit-packing alone. The remaining gap is because:

1. We use uint32 (4 bytes) for delta samples rather than true 32-bit-only mode
2. The delta scheme adds sampling overhead (absolute values every 64 entries)
3. phi_da uses bit_width+1 bits for the -1 sentinel

The achieved compression is substantial and the format is simple, robust, and fully backward-compatible.
