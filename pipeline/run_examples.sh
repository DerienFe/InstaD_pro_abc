#!/bin/bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA="$ROOT/data"
OUT="$ROOT/proteinmpnn_run/results/custom_scores"

python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/6m0j.pdb" \
  --chains "A E" \
  --out-pdb "$DATA/6m0j_ae.pdb" \
  --out-fasta "$DATA/6m0j_ae.fasta" \
  --out-seq "$DATA/6m0j_ae.seq"

python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/7z0x.pdb" \
  --chains "H L R" \
  --out-pdb "$DATA/7z0x_hlr.pdb" \
  --out-fasta "$DATA/7z0x_hlr.fasta" \
  --out-seq "$DATA/7z0x_hlr.seq"

python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_ae.pdb" \
  --design-chains "A E" \
  --fasta "$DATA/6m0j_ae.seq" \
  --design-ranges "A:19-29,324-326,353-355" \
  --num-samples 3 \
  --out-dir "$OUT/6m0j_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_ae.pdb" \
  --design-chains "A" \
  --fasta "$DATA/6m0j_ae.seq" \
  --design-ranges "A:19-29,324-326,353-355" \
  --num-samples 3 \
  --out-dir "$OUT/6m0j_unbound"

python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/7z0x_hlr.pdb" \
  --design-chains "H L R" \
  --fasta "$DATA/7z0x_hlr.seq" \
  --num-samples 1 \
  --out-dir "$OUT/7z0x"
