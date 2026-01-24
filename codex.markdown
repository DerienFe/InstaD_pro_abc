Pipeline blueprint (ProteinMPNN score-only)

Objectives
- Trim PDBs to selected chains and emit chain-order '/'-joined sequences with residue-count validation.
- Score designable residues with ProteinMPNN, capturing per-residue logits/NLL for downstream mutation scans (mask-one-by-one).
- Build antibody numbering maps (Chothia) with CDR annotations to align design positions to structural indices.
- Plot contact maps with dual axes (PDB + Chothia) for antibody complexes.
- Compare bound vs unbound complexes (6m0j_AE vs 6m0j_A) over design ranges 19-29, 324-326, 353-355.

Key scripts (under InstaD_pro_abc)
- preprocessing: pipeline/preprocess_chains.py
- scoring with per-residue logits: pipeline/mpnn_score_only.py
- numbering: pipeline/build_antibody_numbering.py
- plotting: pipeline/plot_contact_map.py
- example driver: pipeline/run_examples.sh
- outputs: proteinmpnn_run/results/custom_scores/

Usage snippets
- Prep 6m0j A/E: python pipeline/preprocess_chains.py --pdb data/6m0j.pdb --chains "A E" --out-pdb data/6m0j_ae.pdb --out-fasta data/6m0j_ae.fasta --out-seq data/6m0j_ae.seq
- Prep 7z0x H/L/R: python pipeline/preprocess_chains.py --pdb data/7z0x.pdb --chains "H L R" --out-pdb data/7z0x_hlr.pdb --out-fasta data/7z0x_hlr.fasta --out-seq data/7z0x_hlr.seq
- Build 7z0x Chothia maps: python pipeline/build_antibody_numbering.py --pdb data/7z0x_hl.pdb --chains "H L" --heavy H --light L --out-csv data/7z0x_chothia_map.csv --out-json data/7z0x_chothia_map.json
- Score 6m0j bound: python pipeline/mpnn_score_only.py --pdb data/6m0j_ae.pdb --design-chains "A E" --fasta data/6m0j_ae.seq --design-ranges "A:19-29,324-326,353-355" --num-samples 3 --ala-scan --all19-scan --out-dir proteinmpnn_run/results/custom_scores/6m0j_bound
- Score 6m0j unbound: python pipeline/mpnn_score_only.py --pdb data/6m0j_ae.pdb --design-chains "A" --fasta data/6m0j_a.seq --design-ranges "A:19-29,324-326,353-355" --num-samples 3 --out-dir proteinmpnn_run/results/custom_scores/6m0j_unbound
- Plot contacts (dual axes): python pipeline/plot_contact_map.py --pdb-a data/7z0x_hl.pdb --chain-a H --pdb-b data/7z0x_hl.pdb --chains-b "H L" --numbering-map data/7z0x_chothia_map.csv --output proteinmpnn_run/results/custom_scores/7z0x_contact.png
- Run all examples: bash pipeline/run_examples.sh

Scoring behavior
- Input sequences: headerless '/'-delimited string for designed chains (ordered as provided after sorting); lengths must match chain lengths. If provided sequences are longer than the observed PDB length, pass --allow-longer-seqs to truncate (useful for antibody numbering gaps).
- Design mask: ranges are 1-indexed inclusive; if omitted, full designed chains are used. Visible chains are fixed automatically.
- Outputs: npz with log_probs [samples,L,21], per_residue_nll [samples,L], native_mean_nll, design_nll, design_mask, global_mask, chain_order, optional alanine and all-19 scan tables; txt with mean design/global NLL plus per-residue and scan summaries.
- Model defaults: vanilla v_48_020, CPU/GPU auto; switch with --ca-only or --use-soluble-model.

Next steps for mutation scans
- For each residue in design ranges, use built-in alanine scan or all-19 scan outputs to assess penalties; compare bound vs unbound npz to prioritize critical positions; align back to Chothia numbering via the generated maps.
