#!/usr/bin/env python
"""Plot residue–residue contact maps between two proteins.

Features:
- CA–CA distance contacts with configurable cutoff.
- Optional Chothia renumbering for the second Prot via ANARCI.
- Crops the map to only contacting regions with a configurable buffer.
"""

import argparse
import os
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


def _residue_key(chain_id: str, residue) -> Tuple[str, int, str]:
    icode = residue.id[2].strip() if residue.id[2].strip() else ""
    return (chain_id, int(residue.id[1]), icode)


def _get_chain(structure, chain_id: str | None):
    chains = list(structure.get_chains())
    if not chains:
        raise ValueError("No chains found in PDB")
    if chain_id is None:
        return chains[0]
    for chain in chains:
        if chain.id == chain_id:
            return chain
    raise ValueError(f"Chain {chain_id} not found; available: {[c.id for c in chains]}")


def _get_chain_required(structure, chain_id: str | None, position_label: str):
    chains = list(structure.get_chains())
    if not chains:
        raise ValueError("No chains found in PDB")
    if chain_id is None:
        if len(chains) < 2:
            raise ValueError(
                f"Need at least two chains to auto-assign {position_label} but only found {[c.id for c in chains]}"
            )
        return chains[0 if position_label.lower().startswith("heavy") else 1]
    for chain in chains:
        if chain.id == chain_id:
            return chain
    raise ValueError(
        f"Chain {chain_id} (for {position_label}) not found; available: {[c.id for c in chains]}"
    )


def _chain_residues_and_sequence(chain) -> Tuple[List, str]:
    residues = []
    seq: List[str] = []
    for res in chain:
        if res.id[0] != " ":
            continue
        if "CA" not in res:
            continue
        try:
            aa = seq1(res.resname)
        except Exception:
            continue
        residues.append(res)
        seq.append(aa)
    return residues, "".join(seq)


def _build_chothia_map(chain) -> Dict[Tuple[str, int, str], str]:
    try:
        from anarci import anarci
    except ImportError as exc:
        raise ImportError(
            "anarci is required for --b-is-ab. Install with: pip install anarci"
        ) from exc

    residues, seq = _chain_residues_and_sequence(chain)
    if not residues:
        raise ValueError("No standard residues with CA atoms found for numbering")

    numbering, _, _ = anarci([(chain.id or "chain", seq)], scheme="chothia")
    if not numbering or numbering[0] is None or numbering[0][0] is None:
        raise ValueError("ANARCI could not number the provided sequence")

    domain_numbering, start_idx, _ = numbering[0][0]
    mapping: Dict[Tuple[str, int, str], str] = {}
    for idx, ((pos, insertion), _aa) in enumerate(domain_numbering):
        seq_idx = start_idx + idx
        if seq_idx >= len(residues):
            break
        label = f"{pos}{insertion.strip()}" if insertion.strip() else str(pos)
        res = residues[seq_idx]
        mapping[_residue_key(chain.id, res)] = label
    return mapping


def _extract_ca(chain, numbering_map: Dict[Tuple[str, int, str], str] | None, prefix: str | None):
    labels: List[str] = []
    coords: List[np.ndarray] = []
    for res in chain:
        if res.id[0] != " ":
            continue
        if "CA" not in res:
            continue
        key = _residue_key(chain.id, res)
        label = numbering_map.get(key) if numbering_map else None
        if label is None:
            resseq = res.id[1]
            icode = res.id[2].strip()
            label = f"{resseq}{icode}" if icode else str(resseq)
        if prefix:
            label = f"{prefix}{label}"
        labels.append(label)
        coords.append(res["CA"].get_coord())
    if not coords:
        raise ValueError(f"No CA atoms found in chain {chain.id}")
    return labels, np.vstack(coords)


def _compute_contacts(coords_a: np.ndarray, coords_b: np.ndarray) -> np.ndarray:
    diff = coords_a[:, None, :] - coords_b[None, :, :]
    return np.linalg.norm(diff, axis=2)


def _select_indices(mask_any: np.ndarray, buffer: int, max_len: int) -> np.ndarray:
    core = np.where(mask_any)[0]
    if core.size == 0:
        return core
    keep = set()
    for idx in core:
        start = max(0, idx - buffer)
        end = min(max_len, idx + buffer + 1)
        for j in range(start, end):
            keep.add(j)
    return np.array(sorted(keep), dtype=int)


def _crop_matrix(dist: np.ndarray, cutoff: float, buffer: int):
    contact_mask = dist <= cutoff
    if not contact_mask.any():
        raise ValueError("No contacts found under the cutoff")
    row_idx = _select_indices(contact_mask.any(axis=1), buffer, dist.shape[0])
    col_idx = _select_indices(contact_mask.any(axis=0), buffer, dist.shape[1])
    return dist[np.ix_(row_idx, col_idx)], (row_idx, col_idx)


def plot_contact_map(
    pdb_a: str,
    pdb_b: str,
    chain_a: str | None,
    chain_b: str | None,
    chain_b_heavy: str | None,
    chain_b_light: str | None,
    cutoff: float,
    buffer: int,
    b_is_ab: bool,
    output: str,
    title: str | None,
):
    parser = PDBParser(QUIET=True)
    structure_a = parser.get_structure("A", pdb_a)
    structure_b = parser.get_structure("B", pdb_b)
    chain_a_obj = _get_chain(structure_a, chain_a)

    if b_is_ab:
        heavy_chain = _get_chain_required(structure_b, chain_b_heavy, "heavy chain")
        light_chain = _get_chain_required(structure_b, chain_b_light, "light chain")
        numbering_map_h = _build_chothia_map(heavy_chain)
        numbering_map_l = _build_chothia_map(light_chain)
        labels_h, coords_h = _extract_ca(heavy_chain, numbering_map_h, "H")
        labels_l, coords_l = _extract_ca(light_chain, numbering_map_l, "L")
        labels_b = labels_h + labels_l
        coords_b = np.vstack([coords_h, coords_l])
    else:
        chain_b_obj = _get_chain(structure_b, chain_b)
        labels_b, coords_b = _extract_ca(chain_b_obj, None, None)

    labels_a, coords_a = _extract_ca(chain_a_obj, None, None)

    distances = _compute_contacts(coords_a, coords_b)
    cropped, (row_idx, col_idx) = _crop_matrix(distances, cutoff, buffer)
    labels_a_crop = [labels_a[i] for i in row_idx]
    labels_b_crop = [labels_b[j] for j in col_idx]

    masked = np.where(cropped <= cutoff, cropped, np.nan)
    cmap = plt.cm.coolwarm
    cmap.set_bad("white")

    fig_w = max(6.0, len(labels_b_crop) * 0.25)
    fig_h = max(5.0, len(labels_a_crop) * 0.25)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(masked, origin="lower", cmap=cmap, aspect="auto")
    cbar = fig.colorbar(im, ax=ax, label="CA-CA distance (Å)")
    cbar.ax.tick_params(labelsize=10)

    ax.set_xticks(np.arange(len(labels_b_crop)))
    ax.set_yticks(np.arange(len(labels_a_crop)))
    ax.set_xticklabels(labels_b_crop, rotation=90)
    ax.set_yticklabels(labels_a_crop)
    ax.set_xlabel("Residue (Prot B)")
    ax.set_ylabel("Residue (Prot A)")
    name_a = os.path.basename(pdb_a)
    name_b = os.path.basename(pdb_b)
    ax.set_title(title or f"Contact map: {name_a} vs {name_b}")
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    print(f"Saved contact map to {output}")


def main(argv: Sequence[str] | None = None):
    ap = argparse.ArgumentParser(description="Plot cropped CA contact map between two PDBs")
    ap.add_argument("pdb_a", help="Path to the first PDB (e.g., pro_a.pdb)")
    ap.add_argument("pdb_b", help="Path to the second PDB (e.g., pro_b.pdb or pro_c.pdb)")
    ap.add_argument("--chain-a", dest="chain_a", default=None, help="Chain ID for Prot A (default: first chain)")
    ap.add_argument("--chain-b", dest="chain_b", default=None, help="Chain ID for Prot B when not antibody (default: first chain)")
    ap.add_argument("--chain-b-heavy", dest="chain_b_heavy", default=None, help="Heavy chain ID for antibody Prot B (default: auto first chain)")
    ap.add_argument("--chain-b-light", dest="chain_b_light", default=None, help="Light chain ID for antibody Prot B (default: auto second chain)")
    ap.add_argument("--cutoff", type=float, default=8.0, help="CA-CA distance cutoff in Å (default: 8.0)")
    ap.add_argument("--buffer", type=int, default=2, help="Buffer (in residues) around contacting regions (default: 10)")
    ap.add_argument("--b-is-ab", action="store_true", help="Treat Prot B as antibody and renumber with Chothia scheme via ANARCI")
    ap.add_argument("--output", default="contact_map.png", help="Output image path (default: contact_map.png)")
    ap.add_argument("--title", default=None, help="Plot title override")
    args = ap.parse_args(argv)

    plot_contact_map(
        pdb_a=args.pdb_a,
        pdb_b=args.pdb_b,
        chain_a=args.chain_a,
        chain_b=args.chain_b,
        chain_b_heavy=args.chain_b_heavy,
        chain_b_light=args.chain_b_light,
        cutoff=args.cutoff,
        buffer=args.buffer,
        b_is_ab=args.b_is_ab,
        output=args.output,
        title=args.title,
    )


if __name__ == "__main__":
    main()
