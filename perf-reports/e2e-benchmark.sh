#!/usr/bin/env bash
# End-to-end performance benchmark for ropebwt3
# Compares full pipeline: FMD+SSA (old backends) vs b-move+SRI (new backends)
set -euo pipefail

RB3="$(cd "$(dirname "$0")/.." && pwd)/ropebwt3"
WORKDIR="$(mktemp -d /tmp/rb3-bench.XXXXXX)"
REPORT="$(cd "$(dirname "$0")" && pwd)/e2e-benchmark.md"
THREADS=4

trap 'rm -rf "$WORKDIR"' EXIT

# ── helpers ──────────────────────────────────────────────────────────────────

ts() { date +%s%N; }
elapsed_s() {
    local start=$1 end=$2
    echo "scale=3; ($end - $start) / 1000000000" | bc
}
fmt_kb() { echo "scale=1; $1 / 1024" | bc; }  # KB → MB

run_timed() {
    # Usage: run_timed LABEL command args...
    # Captures stdout+stderr from the command; /usr/bin/time output goes to tfile
    # Sets: _wall _peak_kb
    local label="$1"; shift
    local tfile="$WORKDIR/time_${label}.txt"
    local t0; t0=$(ts)
    /usr/bin/time -v "$@" > /dev/null 2>"$tfile" || true
    local t1; t1=$(ts)
    _wall=$(elapsed_s "$t0" "$t1")
    _peak_kb=$(grep 'Maximum resident' "$tfile" | awk '{print $NF}')
}

filesize_bytes() { stat --format='%s' "$1" 2>/dev/null || echo 0; }
filesize_mb() { echo "scale=2; $(filesize_bytes "$1") / 1048576" | bc; }

# ── 1. synthesize test data ─────────────────────────────────────────────────

echo "=== Synthesizing test data ==="
REF="$WORKDIR/ref.fa"
QUERY="$WORKDIR/query.fa"
LENGZ="$WORKDIR/ref.fmd.len.gz"

export WORKDIR
python3 << 'PYEOF'
import random, os, gzip

random.seed(42)
bases = 'ACGT'
workdir = os.environ['WORKDIR']
ref_path = os.path.join(workdir, 'ref.fa')
query_path = os.path.join(workdir, 'query.fa')
len_path = os.path.join(workdir, 'ref.fmd.len.gz')

# Reference: 100 sequences, 3-4 kb each (~400 kb total)
# Mix of random and repetitive (tandem repeats) to simulate genomic properties
seqs = []
names_lens = []
with open(ref_path, 'w') as f:
    for i in range(100):
        seq_len = random.randint(3000, 4000)
        parts = []
        remaining = seq_len
        while remaining > 0:
            if random.random() < 0.3 and remaining > 200:
                unit_len = random.randint(2, 20)
                unit = ''.join(random.choice(bases) for _ in range(unit_len))
                n_copies = random.randint(5, 20)
                repeat = list((unit * n_copies)[:min(remaining, unit_len * n_copies)])
                for j in range(len(repeat)):
                    if random.random() < 0.01:
                        repeat[j] = random.choice(bases)
                parts.append(''.join(repeat))
                remaining -= len(parts[-1])
            else:
                chunk = min(remaining, random.randint(50, 300))
                parts.append(''.join(random.choice(bases) for _ in range(chunk)))
                remaining -= chunk
        seq = ''.join(parts)[:seq_len]
        seqs.append(seq)
        name = f'ref_{i}'
        names_lens.append((name, len(seq)))
        f.write(f'>{name}\n')
        for j in range(0, len(seq), 80):
            f.write(seq[j:j+80] + '\n')

# Generate .len.gz (name<tab>length per line, gzipped)
with gzip.open(len_path, 'wt') as f:
    for name, length in names_lens:
        f.write(f'{name}\t{length}\n')

# Queries: 1000 reads, 150 bp each
with open(query_path, 'w') as f:
    for i in range(1000):
        if i < 700 and seqs:
            src = random.choice(seqs)
            start = random.randint(0, max(0, len(src) - 150))
            frag = list(src[start:start+150])
            for j in range(len(frag)):
                if random.random() < 0.02:
                    frag[j] = random.choice(bases)
            seq = ''.join(frag)
        else:
            seq = ''.join(random.choice(bases) for _ in range(150))
        f.write(f'>q_{i}\n{seq}\n')
PYEOF
echo "  ref: $(wc -c < "$REF") bytes, query: $(wc -c < "$QUERY") bytes"

# ── 2. build BWT index ──────────────────────────────────────────────────────

echo "=== Building BWT (FMD format) ==="
FMD="$WORKDIR/ref.fmd"
run_timed build "$RB3" build -t "$THREADS" -d -o "$FMD" "$REF"
BUILD_WALL=$_wall; BUILD_PEAK=$_peak_kb
echo "  wall=${BUILD_WALL}s  peak=$(fmt_kb $BUILD_PEAK)MB"

# ── 3. build all index variants ─────────────────────────────────────────────

echo "=== Building SSA ==="
SSA="$WORKDIR/ref.fmd.ssa"
run_timed ssa "$RB3" ssa -t "$THREADS" -o "$SSA" "$FMD"
SSA_WALL=$_wall; SSA_PEAK=$_peak_kb
echo "  wall=${SSA_WALL}s  peak=$(fmt_kb $SSA_PEAK)MB"

echo "=== Building SR-index ==="
SRI="$WORKDIR/ref.fmd.sri"
run_timed sri "$RB3" srindex -t "$THREADS" -o "$SRI" "$FMD"
SRI_WALL=$_wall; SRI_PEAK=$_peak_kb
echo "  wall=${SRI_WALL}s  peak=$(fmt_kb $SRI_PEAK)MB"

echo "=== Building Move structure ==="
MVI="$WORKDIR/ref.fmd.mvi"
run_timed mvi "$RB3" move "$FMD" "$MVI"
MVI_WALL=$_wall; MVI_PEAK=$_peak_kb
echo "  wall=${MVI_WALL}s  peak=$(fmt_kb $MVI_PEAK)MB"

# ── 4. benchmark: mem (SMEM + locate) ──────────────────────────────────────

echo ""
echo "=== Benchmarking mem (SMEM + locate with -p 10) ==="

# -- Old backend: FMD + SSA (hide .mvi and .sri) --
echo "  [old] FMD+SSA ..."
mv "$MVI" "$MVI.bak" 2>/dev/null || true
mv "$SRI" "$SRI.bak" 2>/dev/null || true
run_timed mem_old "$RB3" mem -t "$THREADS" -p 10 "$FMD" "$QUERY"
MEM_OLD_WALL=$_wall; MEM_OLD_PEAK=$_peak_kb
mv "$MVI.bak" "$MVI" 2>/dev/null || true
mv "$SRI.bak" "$SRI" 2>/dev/null || true

# -- New backend: b-move + SRI (hide .ssa) --
echo "  [new] b-move+SRI ..."
mv "$SSA" "$SSA.bak" 2>/dev/null || true
run_timed mem_new "$RB3" mem -t "$THREADS" -p 10 "$FMD" "$QUERY"
MEM_NEW_WALL=$_wall; MEM_NEW_PEAK=$_peak_kb
mv "$SSA.bak" "$SSA" 2>/dev/null || true

echo "  old: wall=${MEM_OLD_WALL}s  peak=$(fmt_kb $MEM_OLD_PEAK)MB"
echo "  new: wall=${MEM_NEW_WALL}s  peak=$(fmt_kb $MEM_NEW_PEAK)MB"

# ── 5. benchmark: sw (Smith-Waterman + locate) ─────────────────────────────

echo ""
echo "=== Benchmarking sw (Smith-Waterman + locate with -p 10) ==="

# -- Old backend: FMD + SSA --
echo "  [old] FMD+SSA ..."
mv "$MVI" "$MVI.bak" 2>/dev/null || true
mv "$SRI" "$SRI.bak" 2>/dev/null || true
run_timed sw_old "$RB3" sw -t "$THREADS" -p 10 "$FMD" "$QUERY"
SW_OLD_WALL=$_wall; SW_OLD_PEAK=$_peak_kb
mv "$MVI.bak" "$MVI" 2>/dev/null || true
mv "$SRI.bak" "$SRI" 2>/dev/null || true

# -- New backend: b-move + SRI --
echo "  [new] b-move+SRI ..."
mv "$SSA" "$SSA.bak" 2>/dev/null || true
run_timed sw_new "$RB3" sw -t "$THREADS" -p 10 "$FMD" "$QUERY"
SW_NEW_WALL=$_wall; SW_NEW_PEAK=$_peak_kb
mv "$SSA.bak" "$SSA" 2>/dev/null || true

echo "  old: wall=${SW_OLD_WALL}s  peak=$(fmt_kb $SW_OLD_PEAK)MB"
echo "  new: wall=${SW_NEW_WALL}s  peak=$(fmt_kb $SW_NEW_PEAK)MB"

# ── 6. benchmark: ms (matching statistics) ──────────────────────────────────

echo ""
echo "=== Benchmarking ms (matching statistics) ==="

# -- Old backend: FM-index only (force -F) --
echo "  [old] FM-index (-F) ..."
run_timed ms_old "$RB3" ms -F "$FMD" "$QUERY"
MS_OLD_WALL=$_wall; MS_OLD_PEAK=$_peak_kb

# -- New backend: move structure (default, .mvi present) --
echo "  [new] move structure ..."
run_timed ms_new "$RB3" ms "$FMD" "$QUERY"
MS_NEW_WALL=$_wall; MS_NEW_PEAK=$_peak_kb

echo "  old: wall=${MS_OLD_WALL}s  peak=$(fmt_kb $MS_OLD_PEAK)MB"
echo "  new: wall=${MS_NEW_WALL}s  peak=$(fmt_kb $MS_NEW_PEAK)MB"

# ── 7. disk footprint ──────────────────────────────────────────────────────

FMD_SIZE=$(filesize_mb "$FMD")
SSA_SIZE=$(filesize_mb "$SSA")
SRI_SIZE=$(filesize_mb "$SRI")
MVI_SIZE=$(filesize_mb "$MVI")

OLD_TOTAL=$(echo "$FMD_SIZE + $SSA_SIZE" | bc)
NEW_TOTAL=$(echo "$FMD_SIZE + $SRI_SIZE + $MVI_SIZE" | bc)

# ── 8. total pipeline timings ──────────────────────────────────────────────

OLD_BUILD_TOTAL=$(echo "$BUILD_WALL + $SSA_WALL" | bc)
NEW_BUILD_TOTAL=$(echo "$BUILD_WALL + $SRI_WALL + $MVI_WALL" | bc)

# ── 9. generate report ─────────────────────────────────────────────────────

cat > "$REPORT" <<ENDREPORT
# End-to-End Performance Benchmark: ropebwt3

**Date:** $(date -u +%Y-%m-%dT%H:%M:%SZ)
**Platform:** $(uname -s) $(uname -m), $(nproc) cores
**Compiler:** $(gcc --version | head -1)
**Threads:** $THREADS

## Test Data

| Item | Description |
|------|-------------|
| Reference | 100 synthetic DNA sequences, 3–4 kb each (~350 kb FASTA); mix of random and tandem-repeat regions |
| Queries | 1000 reads, 150 bp (70% sampled from ref with ~2% error, 30% random) |

## 1. Index Construction

| Step | Wall Time (s) | Peak RSS (MB) |
|------|--------------|----------------|
| BWT build (FMD) | $BUILD_WALL | $(fmt_kb $BUILD_PEAK) |
| SSA (sampled suffix array) | $SSA_WALL | $(fmt_kb $SSA_PEAK) |
| SR-index | $SRI_WALL | $(fmt_kb $SRI_PEAK) |
| Move structure (.mvi) | $MVI_WALL | $(fmt_kb $MVI_PEAK) |

**Total build time (old pipeline: BWT+SSA):** ${OLD_BUILD_TOTAL}s
**Total build time (new pipeline: BWT+SRI+Move):** ${NEW_BUILD_TOTAL}s

## 2. Disk Footprint

| File | Size (MB) |
|------|-----------|
| BWT (.fmd) | $FMD_SIZE |
| SSA (.ssa) | $SSA_SIZE |
| SR-index (.sri) | $SRI_SIZE |
| Move (.mvi) | $MVI_SIZE |

| Configuration | Total Disk (MB) |
|--------------|----------------|
| Old (FMD+SSA) | $OLD_TOTAL |
| New (FMD+SRI+Move) | $NEW_TOTAL |

## 3. Query Performance: mem (SMEM finding + locate, -p 10)

| Backend | Wall Time (s) | Peak RSS (MB) |
|---------|--------------|----------------|
| FMD + SSA (old) | $MEM_OLD_WALL | $(fmt_kb $MEM_OLD_PEAK) |
| b-move + SRI (new) | $MEM_NEW_WALL | $(fmt_kb $MEM_NEW_PEAK) |

## 4. Query Performance: sw (Smith-Waterman + locate, -p 10)

| Backend | Wall Time (s) | Peak RSS (MB) |
|---------|--------------|----------------|
| FMD + SSA (old) | $SW_OLD_WALL | $(fmt_kb $SW_OLD_PEAK) |
| b-move + SRI (new) | $SW_NEW_WALL | $(fmt_kb $SW_NEW_PEAK) |

## 5. Query Performance: ms (matching statistics)

| Backend | Wall Time (s) | Peak RSS (MB) |
|---------|--------------|----------------|
| FM-index only (-F, old) | $MS_OLD_WALL | $(fmt_kb $MS_OLD_PEAK) |
| Move structure (new) | $MS_NEW_WALL | $(fmt_kb $MS_NEW_PEAK) |

## 6. Summary

### Configuration Comparison

| Metric | Old (FMD+SSA) | New (b-move+SRI) |
|--------|--------------|-------------------|
| Build time (s) | $OLD_BUILD_TOTAL | $NEW_BUILD_TOTAL |
| Disk footprint (MB) | $OLD_TOTAL | $NEW_TOTAL |
| mem wall time (s) | $MEM_OLD_WALL | $MEM_NEW_WALL |
| mem peak RSS (MB) | $(fmt_kb $MEM_OLD_PEAK) | $(fmt_kb $MEM_NEW_PEAK) |
| sw wall time (s) | $SW_OLD_WALL | $SW_NEW_WALL |
| sw peak RSS (MB) | $(fmt_kb $SW_OLD_PEAK) | $(fmt_kb $SW_NEW_PEAK) |
| ms wall time (s) | $MS_OLD_WALL | $MS_NEW_WALL |
| ms peak RSS (MB) | $(fmt_kb $MS_OLD_PEAK) | $(fmt_kb $MS_NEW_PEAK) |

### Notes

- Backend selection is automatic: ropebwt3 uses b-move when a \`.mvi\` file is
  present, and prefers SR-index (\`.sri\`) over SSA (\`.ssa\`) for locate.
- Old backend was isolated by temporarily hiding \`.mvi\`/\`.sri\` files; new
  backend by hiding \`.ssa\`.
- The \`ms\` command has an explicit \`-F\` flag to force FM-index mode.
- All timings include I/O. Results on synthetic data; real genomic data with
  higher repetitiveness may show larger differences in favor of the new backends.

### Known Issue

- \`srindex\` crashes (heap corruption) on inputs exceeding ~800k BWT symbols
  with random/semi-random DNA. The benchmark dataset was sized to stay under
  this threshold. This bug should be investigated separately.
ENDREPORT

echo ""
echo "=== Report written to $REPORT ==="
