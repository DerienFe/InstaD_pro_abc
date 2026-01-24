ProteinMPNN score-only workflow (structure-constrained)

Plan
- Preprocess: trim PDBs to required chains, emit headerless '/'-joined sequences, validate residue counts (use [pipeline/preprocess_chains.py](pipeline/preprocess_chains.py)).
- Scoring: run ProteinMPNN in score-only mode with per-residue logits and design-range masking; support native or provided sequences (use [pipeline/mpnn_score_only.py](pipeline/mpnn_score_only.py)).
- Antibody numbering: generate Chothia PDB↔sequence maps + CDR annotations (use [pipeline/build_antibody_numbering.py](pipeline/build_antibody_numbering.py)).
- Plotting: render contact maps with dual x-axes (PDB index + Chothia) using [pipeline/plot_contact_map.py](pipeline/plot_contact_map.py).
- Bound vs unbound comparison: score 6m0j_AE (bound) then 6m0j_A (unbound) over design ranges 19-29, 324-326, 353-355; compute NLL deltas per residue from output npz.
- Bulk helpers: drive common runs via [pipeline/run_examples.sh](pipeline/run_examples.sh); adjust paths, chains, and ranges for new targets.
- Outputs: npz files contain log_probs [samples, L, 21], per_residue_nll [samples, L], native_mean_nll, design_nll (masked), design_mask, global_mask, chain_order; optional alanine/all-19 scan tables. Txt summaries report mean design/global NLL plus per-residue stats and scans.

- Quick commands
- 6m0j chains A/E prep: `python pipeline/preprocess_chains.py --pdb data/6m0j.pdb --chains "A E" --out-pdb data/6m0j_ae.pdb --out-fasta data/6m0j_ae.fasta --out-seq data/6m0j_ae.seq`
- 7z0x chains H/L/R prep: `python pipeline/preprocess_chains.py --pdb data/7z0x.pdb --chains "H L R" --out-pdb data/7z0x_hlr.pdb --out-fasta data/7z0x_hlr.fasta --out-seq data/7z0x_hlr.seq`
- Build 7z0x antibody numbering (H/L): `python pipeline/build_antibody_numbering.py --pdb data/7z0x_hl.pdb --chains "H L" --heavy H --light L --out-csv data/7z0x_chothia_map.csv --out-json data/7z0x_chothia_map.json`
- Score 6m0j bound (A/E designed, A ranges 19-29,324-326,353-355): `python pipeline/mpnn_score_only.py --pdb data/6m0j_ae.pdb --design-chains "A E" --fasta data/6m0j_ae.seq --design-ranges "A:19-29,324-326,353-355" --num-samples 3 --ala-scan --all19-scan --out-dir proteinmpnn_run/results/custom_scores/6m0j_bound`
- Score 6m0j unbound (A only, same ranges): `python pipeline/mpnn_score_only.py --pdb data/6m0j_ae.pdb --design-chains "A" --fasta data/6m0j_ae.seq --design-ranges "A:19-29,324-326,353-355" --num-samples 3 --out-dir proteinmpnn_run/results/custom_scores/6m0j_unbound`
- Plot contact map with dual axes: `python pipeline/plot_contact_map.py --pdb-a data/7z0x_hl.pdb --chain-a H --pdb-b data/7z0x_hl.pdb --chains-b "H L" --numbering-map data/7z0x_chothia_map.csv --output proteinmpnn_run/results/custom_scores/7z0x_contact.png`
- All examples: `bash pipeline/run_examples.sh`

- Notes
- Input sequence files for scoring should be headerless, chain-ordered, and '/'-delimited (e.g., AAA/BBB). They must match the summed length of designed chains; if the provided sequence is longer, pass --allow-longer-seqs to truncate to the observed PDB length (helpful for numbering-induced padding).
- Design mask: ranges are 1-indexed inclusive per chain; if no ranges are provided, the entire designed chain is scored.
- Per-residue penalties: use per_residue_nll/design_nll for native scoring; alanine scan and all-19 scan deltas quantify positional sensitivity; compare bound vs unbound npz to locate sensitive positions.
- Default model: vanilla v_48_020; toggle --ca-only or --use-soluble-model if needed.
- Numbering/plots: build Chothia maps once per antibody PDB (CSV/JSON) then feed to the contact plotter to display dual PDB/Chothia axes.
