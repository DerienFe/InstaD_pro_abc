#!/bin/bash
set -euo pipefail
ROOT="/home/tiejun/WORK/InstaD_pro_abc"
DATA="$ROOT/data"
OUT="$ROOT/proteinmpnn_run/results/custom_scores"

#S1-ACE2 complex
python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/6m0j.pdb" \
  --chains "A E" \
  --out-pdb "$DATA/6m0j_ae.pdb" \
  --out-fasta "$DATA/6m0j_ae.fasta" \
  --out-seq "$DATA/6m0j_ae.seq"
  
#ACE2 only
python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/6m0j.pdb" \
  --chains "A" \
  --out-pdb "$DATA/6m0j_a.pdb" \
  --out-fasta "$DATA/6m0j_a.fasta" \
  --out-seq "$DATA/6m0j_a.seq"

#S1-antibody complex
python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/7z0x.pdb" \
  --chains "H L R" \
  --out-pdb "$DATA/7z0x_hlr.pdb" \
  --out-fasta "$DATA/7z0x_hlr.fasta" \
  --out-seq "$DATA/7z0x_hlr.seq"

#antibody only
python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/7z0x.pdb" \
  --chains "H L" \
  --out-pdb "$DATA/7z0x_hl.pdb" \
  --out-fasta "$DATA/7z0x_hl.fasta" \
  --out-seq "$DATA/7z0x_hl.seq"

##############################
### MPNN scoring only
##############################

python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_ae.pdb" \
  --design-chains "A E" \
  --fasta "$DATA/6m0j_ae.seq" \
  --design-ranges "A:19-29,324-326,353-355" \
  --num-samples 3 \
  --out-dir "$OUT/6m0j_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_a.pdb" \
  --design-chains "A" \
  --fasta "$DATA/6m0j_a.seq" \
  --design-ranges "A:19-29,324-326,353-355" \
  --num-samples 3 \
  --out-dir "$OUT/6m0j_unbound"


#antibody case
python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/7z0x_hlr.pdb" \
  --design-chains "H L R" \
  --fasta "$DATA/7z0x_hlr.seq" \
   --design-ranges "H:30-34,50-58,97-100"\
   --allow-longer-seqs \
  --num-samples 100 \
  --out-dir "$OUT/7z0x_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
    --pdb "$DATA/7z0x_hl.pdb" \
    --design-chains "H L" \
    --fasta "$DATA/7z0x_hl.seq" \
    --design-ranges "H:30-34,50-58,97-100"\
    --allow-longer-seqs \
    --num-samples 100 \
    --out-dir "$OUT/7z0x_unbound"