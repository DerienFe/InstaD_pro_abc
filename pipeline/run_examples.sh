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

#S1 only (ACE2-S1 complex)
python "$ROOT/pipeline/preprocess_chains.py" \
  --pdb "$DATA/6m0j.pdb" \
  --chains "E" \
  --out-pdb "$DATA/6m0j_e.pdb" \
  --out-fasta "$DATA/6m0j_e.fasta" \
  --out-seq "$DATA/6m0j_e.seq"


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

#S1 only (S1-antibody complex)
python "$ROOT/pipeline/preprocess_chains.py" \
     --pdb "$DATA/7z0x.pdb" \
    --chains "R" \
    --out-pdb "$DATA/7z0x_r.pdb" \
    --out-fasta "$DATA/7z0x_r.fasta" \
    --out-seq "$DATA/7z0x_r.seq"

##############################
######## MPNN scoring ########
##############################

#note:
# below are 4 sets of MPNN scoring runs:
# 1. ACE2-S1 complex vs ACE2 only
# 2. ACE2-S1 complex vs S1 only
# 3. Antibody-S1 complex vs Antibody only
# 4. Antibody-S1 complex vs S1 only
# the idea is similar to MM-GBSA binding energy estimation:
# structure constrainted design of the binding partner in complex vs unbound state
# the difference in design NLL is expected to reflect the binding contribution of the designed region
##################################

#first we run ACE2-S1 complex and ACE2 only
#ACE2 side
python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_ae.pdb" \
  --design-chains "A E" \
  --fasta "$DATA/6m0j_ae.seq" \
  --design-ranges "A:19-29,324-326,353-355" \
  --num-samples 100 \
  --out-dir "$OUT/6m0j_ace_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_a.pdb" \
  --design-chains "A" \
  --fasta "$DATA/6m0j_a.seq" \
  --design-ranges "A:19-29,324-326,353-355" \
  --num-samples 100 \
  --out-dir "$OUT/6m0j_ace_unbound"

#S1 side
python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/6m0j_ae.pdb" \
  --design-chains "A E" \
  --fasta "$DATA/6m0j_ae.seq" \
    --design-ranges "E:412-420,452-457" \
    --num-samples 100 \
    --out-dir "$OUT/6m0j_s1_only_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
    --pdb "$DATA/6m0j_e.pdb" \
    --design-chains "E" \
    --fasta "$DATA/6m0j_e.seq" \
    --design-ranges "E:412-420,452-457" \
    --num-samples 100 \
    --out-dir "$OUT/6m0j_s1_only_unbound"

##################################
#antibody case
python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/7z0x_hlr.pdb" \
  --design-chains "H L R" \
  --fasta "$DATA/7z0x_hlr.seq" \
   --design-ranges "H:26-33,51-60,98-116;L:24-37,92-101"\
   --allow-longer-seqs \
  --num-samples 100 \
  --out-dir "$OUT/7z0x_ab_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
    --pdb "$DATA/7z0x_hl.pdb" \
    --design-chains "H L" \
    --fasta "$DATA/7z0x_hl.seq" \
    --design-ranges "H:26-33,51-60,98-116;L:24-37,92-101"\
    --allow-longer-seqs \
    --num-samples 100 \
    --out-dir "$OUT/7z0x_ab_unbound"

#S1 side
#note we have to redesignate the desinable region to S1 range
python "$ROOT/pipeline/mpnn_score_only.py" \
  --pdb "$DATA/7z0x_hlr.pdb" \
    --design-chains "H L R" \
    --fasta "$DATA/7z0x_hlr.seq" \
    --design-ranges "R:412-420,452-457" \
    --num-samples 100 \
    --out-dir "$OUT/7z0x_s1_only_bound"

python "$ROOT/pipeline/mpnn_score_only.py" \
    --pdb "$DATA/7z0x_r.pdb" \
    --design-chains "R" \
    --fasta "$DATA/7z0x_r.seq" \
    --design-ranges "R:412-420,452-457" \
    --num-samples 100 \
    --out-dir "$OUT/7z0x_s1_only_unbound"

##############################
### Delta analysis
##############################

# 6m0j bound vs unbound
#effect on ACE2:
python "$ROOT/pipeline/analyze_deltas.py" \
  --bound-csv "$OUT/6m0j_ace_bound/6m0j_ae_design_nll.csv" \
  --unbound-csv "$OUT/6m0j_ace_unbound/6m0j_a_design_nll.csv" \
  --label 6m0j \
  --map-json "$DATA/6m0j_ae_map.json" \
  --out-dir "$ROOT/proteinmpnn_run/results/analyzed"
#effect on S1:
python "$ROOT/pipeline/analyze_deltas.py" \
  --bound-csv "$OUT/6m0j_s1_only_bound/6m0j_ae_design_nll.csv" \
  --unbound-csv "$OUT/6m0j_s1_only_unbound/6m0j_e_design_nll.csv" \
  --label 6m0j_s1 \
  --map-json "$DATA/6m0j_ae_map.json" \
  --out-dir "$ROOT/proteinmpnn_run/results/analyzed"


# 7z0x bound vs unbound (with Chothia annotation) 
#effect on antibody:
python "$ROOT/pipeline/analyze_deltas.py" \
  --bound-csv "$OUT/7z0x_ab_bound/7z0x_hlr_design_nll.csv" \
  --unbound-csv "$OUT/7z0x_ab_unbound/7z0x_hl_design_nll.csv" \
  --label 7z0x \
  --chothia-map "$ROOT/data/7z0x_chothia_map.json" \
  --out-dir "$ROOT/proteinmpnn_run/results/analyzed"
#effect on S1:
python "$ROOT/pipeline/analyze_deltas.py" \
  --bound-csv "$OUT/7z0x_s1_only_bound/7z0x_hlr_design_nll.csv" \
  --unbound-csv "$OUT/7z0x_s1_only_unbound/7z0x_r_design_nll.csv" \
  --label 7z0x_s1 \
  --map-json "$DATA/7z0x_hlr_map.json" \
  --out-dir "$ROOT/proteinmpnn_run/results/analyzed"