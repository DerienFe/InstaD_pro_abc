#!/usr/bin/env python3
"""
Score-only driver for ProteinMPNN with per-residue log-probs.
- Uses ProteinMPNN directly to emit residue-wise NLL/logits for masked chains.
- Accepts headerless '/'-joined design sequences or FASTA with the same content.
"""
import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import PDB

import numpy as np
import torch
import torch.nn.functional as F


ALPHABET = 'ACDEFGHIKLMNPQRSTVWYX'
ALPHABET_DICT = {aa: idx for idx, aa in enumerate(ALPHABET)}


def parse_design_ranges(raw: str | None) -> Dict[str, List[Tuple[int, int]]]:
    if not raw:
        return {}

    # Accept semicolon-separated chain blocks; if missing, split on repeated "<chain>:" patterns.
    matches = list(re.finditer(r"\s*([A-Za-z0-9]+)\s*:", raw))
    if not matches:
        raise ValueError(f"Invalid design ranges string '{raw}'. Expected format like 'H:10-20,30-35;L:5-9'.")

    blocks: List[str] = []
    for idx, match in enumerate(matches):
        start = match.start()
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(raw)
        block = raw[start:end].strip().rstrip(';,')
        if block:
            blocks.append(block)

    out: Dict[str, List[Tuple[int, int]]] = {}
    for block in blocks:
        chain_part, ranges_part = block.split(':', 1)
        chain = chain_part.strip()
        if not chain:
            raise ValueError(f"Empty chain id in design range block '{block}'.")

        ranges_text = ranges_part.strip()
        if not ranges_text:
            raise ValueError(f"Empty range list for chain '{chain}'.")

        spans: List[Tuple[int, int]] = []
        for rng in ranges_text.split(','):
            rng = rng.strip()
            if not rng:
                continue
            if '-' in rng:
                a_str, b_str = rng.split('-', 1)
                start, end = int(a_str), int(b_str)
            else:
                start = end = int(rng)
            if start > end:
                raise ValueError(f"Range start greater than end for chain '{chain}': '{rng}'.")
            spans.append((start, end))

        if not spans:
            raise ValueError(f"No valid ranges for chain '{chain}' in '{block}'.")
        out[chain] = spans

    return out


def clean_chain_list(raw: str) -> List[str]:
    parts = raw.replace(',', ' ').split()
    uniq = []
    for p in parts:
        if p not in uniq:
            uniq.append(p)
    return sorted(uniq)


def read_design_sequence(path: Path, chain_ids: List[str], chain_lengths: Dict[str, int], allow_longer: bool) -> List[str]:
    text = path.read_text().strip().splitlines()
    text = [line for line in text if not line.startswith('>')]
    seq_joined = ''.join(text).replace(' ', '')
    seq_parts = seq_joined.split('/') if '/' in seq_joined else [seq_joined]
    if len(seq_parts) != len(chain_ids):
        raise ValueError(f"Expected {len(chain_ids)} chain sequences for {chain_ids}, got {len(seq_parts)} in {path.name}")
    cleaned: List[str] = []
    for cid, seq in zip(chain_ids, seq_parts):
        expected = chain_lengths[cid]
        if len(seq) != expected:
            if allow_longer and len(seq) > expected:
                print(f"[warn] Truncating chain {cid} sequence from {len(seq)} to expected {expected} based on PDB length", file=sys.stderr)
                seq = seq[:expected]
            elif len(seq) < expected:
                pad = expected - len(seq)
                print(f"[warn] Padding chain {cid} sequence from {len(seq)} to expected {expected} with {pad} 'X' residues (PDB has internal gaps)", file=sys.stderr)
                seq = seq + ('X' * pad)
            else:
                head = seq[:20]
                tail = seq[-20:] if len(seq) > 20 else seq
                raise ValueError(
                    f"Length mismatch for chain {cid} in {path.name}: expected {expected}, got {len(seq)}. "
                    f"First20={head} | Last20={tail}"
                )
        cleaned.append(seq)
    return cleaned


def build_chain_spans(pdb_dict: dict, design_chains: List[str], visible_chains: List[str]) -> Tuple[Dict[str, Tuple[int, int]], List[str]]:
    chain_order = design_chains + visible_chains
    spans: Dict[str, Tuple[int, int]] = {}
    cursor = 0
    for cid in chain_order:
        seq = pdb_dict[f'seq_chain_{cid}']
        spans[cid] = (cursor, cursor + len(seq))
        cursor += len(seq)
    return spans, chain_order


def load_mpnn(mpnn_root: Path, ca_only: bool, model_name: str, use_soluble: bool) -> torch.nn.Module:
    sys.path.insert(0, str(mpnn_root))
    sys.path.insert(0, str(mpnn_root.resolve()))
    from protein_mpnn_utils import ProteinMPNN  # type: ignore

    weight_dir = mpnn_root / ('ca_model_weights' if ca_only else ('soluble_model_weights' if use_soluble else 'vanilla_model_weights'))
    checkpoint = torch.load(weight_dir / f"{model_name}.pt", map_location='cpu')
    model = ProteinMPNN(
        ca_only=ca_only,
        num_letters=21,
        node_features=128,
        edge_features=128,
        hidden_dim=128,
        num_encoder_layers=3,
        num_decoder_layers=3,
        augment_eps=0.0,
        k_neighbors=checkpoint['num_edges'],
    )
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    return model


def position_lookup(spans: Dict[str, Tuple[int, int]], pdb_indices: Dict[str, List[Tuple[int, str]]]) -> List[Dict[str, object]]:
    positions: List[Dict[str, object]] = []
    for chain, (start, end) in spans.items():
        pdb_list = pdb_indices.get(chain, [])
        for offset in range(start, end):
            local_idx = offset - start
            pdb_idx, icode = pdb_list[local_idx] if local_idx < len(pdb_list) else (local_idx + 1, '')
            positions.append({'chain': chain, 'global_idx': offset, 'chain_idx': local_idx + 1, 'pdb_idx': pdb_idx, 'icode': icode})
    positions = sorted(positions, key=lambda x: x['global_idx'])
    return positions


def pdb_index_map(pdb_path: Path, chain_ids: List[str]) -> Dict[str, List[Tuple[int, str]]]:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, pdb_path)
    mapping: Dict[str, List[Tuple[int, str]]] = {}
    for cid in chain_ids:
        chain = structure[0][cid]
        entries: List[Tuple[int, str]] = []
        for res in chain:
            if res.id[0] != ' ':
                continue
            if 'CA' not in res:
                continue
            resseq = int(res.id[1])
            icode = res.id[2].strip() or ''
            entries.append((resseq, icode))
        mapping[cid] = entries
    return mapping


def tensors_from_pdb(mpnn_root: Path, pdb_path: Path, design_chains: List[str], ca_only: bool, design_ranges: Dict[str, List[Tuple[int, int]]]):
    sys.path.insert(0, str(mpnn_root))
    from protein_mpnn_utils import parse_PDB, StructureDatasetPDB, tied_featurize  # type: ignore

    pdb_dict_list = parse_PDB(str(pdb_path), ca_only=ca_only)
    if not pdb_dict_list:
        raise ValueError(f"Failed to parse {pdb_path}")
    pdb_dict = pdb_dict_list[0]
    all_chains = sorted([k.split('_')[-1] for k in pdb_dict if k.startswith('seq_chain_')])
    visible = [c for c in all_chains if c not in design_chains]
    chain_id_dict = {pdb_dict['name']: (design_chains, visible)}
    dataset = StructureDatasetPDB(pdb_dict_list, truncate=None, max_length=200000)
    protein = dataset[0]
    X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list, visible_list_list, masked_list_list, masked_chain_length_list_list, chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask, tied_pos_list_of_lists_list, pssm_coef_all, pssm_bias_all, pssm_log_odds_all, bias_by_res_all, tied_beta = tied_featurize(
        [protein],
        torch.device('cpu'),
        chain_id_dict,
        fixed_position_dict=None,
        omit_AA_dict=None,
        tied_positions_dict=None,
        pssm_dict=None,
        bias_by_res_dict=None,
        ca_only=ca_only,
    )
    spans, chain_order = build_chain_spans(pdb_dict, design_chains, visible)
    pdb_indices = pdb_index_map(pdb_path, chain_order)
    for cid in chain_order:
        expected = spans[cid][1] - spans[cid][0]
        observed = len(pdb_indices.get(cid, []))
        if expected != observed:
            print(f"[warn] chain {cid}: length {expected} vs pdb-indexed residues {observed}; downstream masks use pdb indices where available", file=sys.stderr)
    design_mask = torch.zeros_like(mask)
    for cid in design_chains:
        start, end = spans[cid]
        length = end - start
        if cid in design_ranges:
            local_mask = torch.zeros(length, device=design_mask.device)
            pdb_list = pdb_indices.get(cid, [])
            for idx in range(length):
                pdb_idx = pdb_list[idx][0] if idx < len(pdb_list) else (idx + 1)
                for lo, hi in design_ranges[cid]:
                    if lo <= pdb_idx <= hi:
                        local_mask[idx] = 1.0
                        break
        else:
            local_mask = torch.ones(length, device=design_mask.device)
        design_mask[:, start:end] = local_mask
    base_mask = mask * chain_M * chain_M_pos
    design_mask = design_mask * base_mask
    positions = position_lookup(spans, pdb_indices)
    return pdb_dict, X, S, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, design_mask, spans, chain_order, positions


def override_sequence(S: torch.Tensor, spans: Dict[str, Tuple[int, int]], design_chains: List[str], seq_parts: List[str]) -> torch.Tensor:
    S_out = S.clone()
    for cid, seq in zip(design_chains, seq_parts):
        start, end = spans[cid]
        encoded = torch.tensor([ALPHABET_DICT[aa] for aa in seq], dtype=torch.long)
        if len(encoded) != (end - start):
            raise ValueError(f"Chain {cid} encoded length mismatch: got {len(encoded)}, expected {end-start}")
        S_out[:, start:end] = encoded
    return S_out


def forward_log_probs(model: torch.nn.Module, X: torch.Tensor, S: torch.Tensor, mask: torch.Tensor, chain_M: torch.Tensor, chain_M_pos: torch.Tensor, residue_idx: torch.Tensor, chain_encoding_all: torch.Tensor, noise: torch.Tensor | None = None) -> torch.Tensor:
    device = next(model.parameters()).device
    X = X.to(device)
    S = S.to(device)
    mask = mask.to(device)
    chain_M = chain_M.to(device)
    chain_M_pos = chain_M_pos.to(device)
    residue_idx = residue_idx.to(device)
    chain_encoding_all = chain_encoding_all.to(device)
    if noise is None:
        noise = torch.zeros_like(chain_M, device=device)
    else:
        noise = noise.to(device)
    log_probs = model(X, S, mask, chain_M * chain_M_pos, residue_idx, chain_encoding_all, noise)
    return log_probs


def score_sequences(model: torch.nn.Module, X: torch.Tensor, S: torch.Tensor, mask: torch.Tensor, chain_M: torch.Tensor, chain_M_pos: torch.Tensor, residue_idx: torch.Tensor, chain_encoding_all: torch.Tensor, design_mask: torch.Tensor, num_samples: int) -> dict:
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    X = X.to(device)
    S = S.to(device)
    mask = mask.to(device)
    chain_M = chain_M.to(device)
    chain_M_pos = chain_M_pos.to(device)
    residue_idx = residue_idx.to(device)
    chain_encoding_all = chain_encoding_all.to(device)
    design_mask = design_mask.to(device)

    per_sample_log_probs = []
    per_sample_nll = []
    per_sample_design = []
    per_sample_global = []

    with torch.no_grad():
        for _ in range(num_samples):
            noise = torch.randn(chain_M.shape, device=device)
            log_probs = forward_log_probs(model, X, S, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, noise)
            nll = F.nll_loss(log_probs[0], S[0], reduction='none')
            design_mean = (nll * design_mask[0]).sum() / design_mask[0].sum()
            global_mean = (nll * mask[0]).sum() / mask[0].sum()
            per_sample_log_probs.append(log_probs.cpu().numpy())
            per_sample_nll.append(nll.cpu().numpy())
            per_sample_design.append(float(design_mean.cpu().numpy()))
            per_sample_global.append(float(global_mean.cpu().numpy()))

    return {
        'log_probs': np.stack(per_sample_log_probs, axis=0),
        'per_residue_nll': np.stack(per_sample_nll, axis=0),
        'design_mean': np.array(per_sample_design, dtype=np.float32),
        'global_mean': np.array(per_sample_global, dtype=np.float32),
    }


def run_ala_scan(model: torch.nn.Module, X: torch.Tensor, S_native: torch.Tensor, mask: torch.Tensor, chain_M: torch.Tensor, chain_M_pos: torch.Tensor, residue_idx: torch.Tensor, chain_encoding_all: torch.Tensor, design_mask: torch.Tensor, native_mean_nll: np.ndarray, positions: List[Dict[str, object]]) -> List[Dict[str, object]]:
    device = next(model.parameters()).device
    results: List[Dict[str, object]] = []
    ala_idx = ALPHABET_DICT['A']
    with torch.no_grad():
        for pos_info in positions:
            gidx = pos_info['global_idx']
            if design_mask[0, gidx] <= 0:
                continue
            S_mut = S_native.clone()
            S_mut[0, gidx] = ala_idx
            log_probs_mut = forward_log_probs(model, X, S_mut, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, noise=None)
            nll_mut = F.nll_loss(log_probs_mut[0], S_mut[0], reduction='none')
            native_nll = float(native_mean_nll[gidx])
            ala_nll = float(nll_mut[gidx].cpu().numpy())
            results.append({
                'chain': pos_info['chain'],
                'chain_idx': pos_info['chain_idx'],
                'global_idx': gidx,
                'native_aa': ALPHABET[int(S_native[0, gidx].cpu().numpy())],
                'native_nll': native_nll,
                'ala_nll': ala_nll,
                'delta': ala_nll - native_nll,
            })
    return results


def run_all19_scan(model: torch.nn.Module, X: torch.Tensor, S_native: torch.Tensor, mask: torch.Tensor, chain_M: torch.Tensor, chain_M_pos: torch.Tensor, residue_idx: torch.Tensor, chain_encoding_all: torch.Tensor, design_mask: torch.Tensor, positions: List[Dict[str, object]]) -> List[Dict[str, object]]:
    device = next(model.parameters()).device
    results: List[Dict[str, object]] = []
    with torch.no_grad():
        log_probs = forward_log_probs(model, X, S_native, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, noise=None)
        log_probs_cpu = log_probs[0].cpu().numpy()  # [L, 21]
        for pos_info in positions:
            gidx = pos_info['global_idx']
            if design_mask[0, gidx] <= 0:
                continue
            native_idx = int(S_native[0, gidx].cpu().numpy())
            if native_idx >= 20:
                continue  # skip unknown residue
            nll_vec = -log_probs_cpu[gidx, :20]  # exclude X
            native_nll = float(nll_vec[native_idx])
            non_native = [i for i in range(20) if i != native_idx]
            mean_non_native = float(nll_vec[non_native].mean())
            results.append({
                'chain': pos_info['chain'],
                'chain_idx': pos_info['chain_idx'],
                'global_idx': gidx,
                'native_aa': ALPHABET[native_idx],
                'native_nll': native_nll,
                'non_native_mean_nll': mean_non_native,
                'delta': mean_non_native - native_nll,
            })
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="ProteinMPNN score-only with per-residue logits")
    parser.add_argument('--pdb', required=True, type=Path)
    parser.add_argument('--design-chains', required=True, help="Chains to score/design (space/comma separated)")
    parser.add_argument('--fasta', type=Path, default=None, help="Headerless '/'-joined sequences for design chains")
    parser.add_argument('--design-ranges', type=str, default=None, help="Semicolon-separated chain blocks; commas between ranges. e.g. H:26-33,51-57,98-116;L:24-37,92-101")
    parser.add_argument('--out-dir', required=True, type=Path)
    parser.add_argument('--num-samples', type=int, default=1000)
    parser.add_argument('--model-name', type=str, default='v_48_020')
    parser.add_argument('--ca-only', action='store_true')
    parser.add_argument('--use-soluble-model', action='store_true')
    parser.add_argument('--mpnn-root', type=Path, default=Path(__file__).resolve().parents[2] / 'ProteinMPNN')
    parser.add_argument('--ala-scan', action='store_true', help='Run per-position alanine scan over designable residues')
    parser.add_argument('--all19-scan', action='store_true', help='Run per-position all-19 average penalty over designable residues')
    parser.add_argument('--allow-longer-seqs', action='store_true', help='Allow provided sequences longer than the PDB length by truncating to observed length')
    args = parser.parse_args()

    design_chains = clean_chain_list(args.design_chains)
    design_ranges = parse_design_ranges(args.design_ranges)

    model = load_mpnn(args.mpnn_root, args.ca_only, args.model_name, args.use_soluble_model)
    pdb_dict, X, S_native, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, design_mask, spans, chain_order, positions = tensors_from_pdb(
        args.mpnn_root,
        args.pdb,
        design_chains,
        args.ca_only,
        design_ranges,
    )

    chain_lengths = {cid: spans[cid][1] - spans[cid][0] for cid in design_chains}
    if args.fasta:
        seq_parts = read_design_sequence(args.fasta, design_chains, chain_lengths, args.allow_longer_seqs)
        S_used = override_sequence(S_native, spans, design_chains, seq_parts)
        seq_override = '/'.join(seq_parts)
        print(f"[info] Using provided design sequence from {args.fasta}: {seq_override}")
    else:
        S_used = S_native
        seq_override = None
        print("[info] Using native sequence from PDB for designed chains")

    scores = score_sequences(model, X, S_used, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, design_mask, args.num_samples)

    native_mean_nll = scores['per_residue_nll'].mean(axis=0)
    design_positions = [p for p in positions if design_mask[0, p['global_idx']] > 0]
    design_nll = native_mean_nll * design_mask[0].cpu().numpy()

    ala_scan = run_ala_scan(model, X, S_used, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, design_mask, native_mean_nll, design_positions) if args.ala_scan else []
    all19_scan = run_all19_scan(model, X, S_used, mask, chain_M, chain_M_pos, residue_idx, chain_encoding_all, design_mask, design_positions) if args.all19_scan else []

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.pdb.stem
    out_npz = args.out_dir / f"{stem}_score_only.npz"
    np.savez(
        out_npz,
        log_probs=scores['log_probs'],
        per_residue_nll=scores['per_residue_nll'],
        design_mean=scores['design_mean'],
        global_mean=scores['global_mean'],
        design_mask=design_mask.cpu().numpy(),
        global_mask=mask.cpu().numpy(),
        chain_order=chain_order,
        design_chains=design_chains,
        sequence_override=seq_override,
        native_mean_nll=native_mean_nll,
        design_nll=design_nll,
        positions=positions,
        ala_scan=ala_scan,
        all19_scan=all19_scan,
    )

    design_vals = design_nll[design_nll > 0]
    design_stats = f"design_residue_NLL mean={design_vals.mean():.4f}, median={np.median(design_vals):.4f}, min={design_vals.min():.4f}, max={design_vals.max():.4f}" if design_vals.size else "design_residue_NLL n/a"

    # Save design-only per-residue NLL to CSV for easier inspection
    design_rows = ["chain,pdb_idx,chain_idx,global_idx,aa,nll"]
    for pos in design_positions:
        gidx = pos['global_idx']
        nll_val = design_nll[gidx]
        aa_idx = int(S_used[0, gidx].cpu().numpy())
        aa = ALPHABET[aa_idx] if aa_idx < len(ALPHABET) else 'X'
        design_rows.append(f"{pos['chain']},{pos.get('pdb_idx','')},{pos['chain_idx']},{gidx},{aa},{nll_val:.6f}")
    (args.out_dir / f"{stem}_design_nll.csv").write_text('\n'.join(design_rows) + '\n')

    out_txt = args.out_dir / f"{stem}_summary.txt"
    lines = [
        f"pdb={args.pdb}",
        f"chains_design={design_chains}",
        f"design_ranges={design_ranges if design_ranges else 'full masked chains'}",
        f"num_samples={args.num_samples}",
        f"design_mean_NLL={scores['design_mean'].mean():.4f} +- {scores['design_mean'].std():.4f}",
        f"global_mean_NLL={scores['global_mean'].mean():.4f} +- {scores['global_mean'].std():.4f}",
        design_stats,
        f"design_nll_values={[float(x) for x in design_nll.tolist()]}",
    ]
    if seq_override:
        lines.append(f"sequence_override={seq_override}")
    if ala_scan:
        lines.append("ala_scan (chain,pos,native,ala_nll,delta):")
        for entry in ala_scan:
            lines.append(f"  {entry['chain']}:{entry['chain_idx']} {entry['native_aa']} -> A | ala_nll={entry['ala_nll']:.4f} delta={entry['delta']:.4f}")
    if all19_scan:
        lines.append("all19_scan (chain,pos,native,mean_non_native_nll,delta):")
        for entry in all19_scan:
            lines.append(f"  {entry['chain']}:{entry['chain_idx']} {entry['native_aa']} | mean_non_native={entry['non_native_mean_nll']:.4f} delta={entry['delta']:.4f}")

    out_txt.write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
