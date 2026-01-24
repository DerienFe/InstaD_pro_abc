Pipeline blueprint (ProteinMPNN score-only)

Objectives
- Trim PDBs to selected chains and emit chain-order '/'-joined sequences with residue-count validation.
- Score designable residues with ProteinMPNN, capturing per-residue logits/NLL for downstream mutation scans (mask-one-by-one).
- Compare bound vs unbound complexes (6m0j_AE vs 6m0j_A) over design ranges 19-29, 324-326, 353-355.

Key scripts (under InstaD_pro_abc)
- preprocessing: pipeline/preprocess_chains.py
- scoring with per-residue logits: pipeline/mpnn_score_only.py
- example driver: pipeline/run_examples.sh
- outputs: proteinmpnn_run/results/custom_scores/

Usage snippets
- Prep 6m0j A/E: python pipeline/preprocess_chains.py --pdb data/6m0j.pdb --chains "A E" --out-pdb data/6m0j_ae.pdb --out-fasta data/6m0j_ae.fasta --out-seq data/6m0j_ae.seq
- Prep 7z0x H/L/R: python pipeline/preprocess_chains.py --pdb data/7z0x.pdb --chains "H L R" --out-pdb data/7z0x_hlr.pdb --out-fasta data/7z0x_hlr.fasta --out-seq data/7z0x_hlr.seq
- Score 6m0j bound: python pipeline/mpnn_score_only.py --pdb data/6m0j_ae.pdb --design-chains "A E" --fasta data/6m0j_ae.seq --design-ranges "A:19-29,324-326,353-355" --num-samples 3 --ala-scan --all19-scan --out-dir proteinmpnn_run/results/custom_scores/6m0j_bound
- Score 6m0j unbound: python pipeline/mpnn_score_only.py --pdb data/6m0j_ae.pdb --design-chains "A" --fasta data/6m0j_ae.seq --design-ranges "A:19-29,324-326,353-355" --num-samples 3 --out-dir proteinmpnn_run/results/custom_scores/6m0j_unbound
- Run all examples: bash pipeline/run_examples.sh

Scoring behavior
- Input sequences: headerless '/'-delimited string for designed chains (ordered as provided after sorting); lengths must match chain lengths.
- Design mask: ranges are 1-indexed inclusive; if omitted, full designed chains are used. Visible chains are fixed automatically.
- Outputs: npz with log_probs [samples,L,21], per_residue_nll [samples,L], native_mean_nll, design_nll, design_mask, global_mask, chain_order, optional alanine and all-19 scan tables; txt with mean design/global NLL plus per-residue and scan summaries.
- Model defaults: vanilla v_48_020, CPU/GPU auto; switch with --ca-only or --use-soluble-model.

Next steps for mutation scans
- For each residue in design ranges, use built-in alanine scan or all-19 scan outputs to assess penalties; compare bound vs unbound npz to prioritize critical positions.
