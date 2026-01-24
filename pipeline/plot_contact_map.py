#!/usr/bin/env python3
"""Contact map plotter with dual axes (PDB index + Chothia labels).

- Computes CA–CA distance contacts between chain A and chain(s) B from PDB files.
- Accepts numbering maps (CSV/JSON) emitted by build_antibody_numbering.py for Chothia labels.
- Crops to contacting regions (with buffer) and renders dual x-axes: bottom = PDB/chain index, top = Chothia if available.
"""
import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

NumberingMap = Dict[str, List[Dict[str, object]]]


def _clean_list(raw: str) -> List[str]:
    parts = raw.replace(",", " ").split()
    uniq: List[str] = []
    for p in parts:
        if p and p not in uniq:
            uniq.append(p)
    return uniq


def _residue_key(res) -> Tuple[int, str]:
    return int(res.id[1]), res.id[2].strip() or ""


def _iter_ca(chain):
    coords: List[np.ndarray] = []
    labels_chain: List[Tuple[int, str, str]] = []  # pdb_idx, icode, aa
    for res in chain:
        if res.id[0] != " ":
            continue
        if "CA" not in res:
            continue
        try:
            aa = seq1(res.resname)
        except Exception:
            continue
        coords.append(res["CA"].get_coord())
        labels_chain.append((int(res.id[1]), res.id[2].strip() or "", aa))
    if not coords:
        raise ValueError(f"No CA atoms found in chain {chain.id}")
    return np.vstack(coords), labels_chain


def _load_numbering_map(path: Path) -> NumberingMap:
    if not path:
        return {}
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text())
        chains = data.get("map", {}) if isinstance(data, dict) else {}
        parsed: NumberingMap = {}
        for cid, payload in chains.items():
            rows = payload.get("rows", []) if isinstance(payload, dict) else []
            parsed[cid] = rows
        return parsed
    parsed: NumberingMap = {}
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        return parsed
    header = lines[0].split(',')
    col_idx = {name: i for i, name in enumerate(header)}
    for line in lines[1:]:
        cols = line.split(',')
        cid = cols[col_idx.get("chain", 0)]
        parsed.setdefault(cid, []).append(
            {
                "chain": cid,
                "chain_idx": int(cols[col_idx.get("chain_idx", 0)]),
                "pdb_idx": int(cols[col_idx.get("pdb_idx", 0)]),
                "icode": cols[col_idx.get("icode", 0)],
                "chothia": cols[col_idx.get("chothia", 0)],
                "aa": cols[col_idx.get("aa", 0)],
                "cdr": cols[col_idx.get("cdr", 0)] if "cdr" in col_idx else "",
            }
        )
    return parsed


def _labels_for_chain(chain_id: str, ca_labels: List[Tuple[int, str, str]], number_map: NumberingMap):
    pdb_labels = []
    chothia_labels = []
    rows = number_map.get(chain_id, [])
    rows_by_idx = {int(r.get("chain_idx", i + 1)): r for i, r in enumerate(rows)}
    for idx, (pdb_idx, icode, _aa) in enumerate(ca_labels, start=1):
        pdb_label = f"{chain_id}{pdb_idx}{icode}" if icode else f"{chain_id}{pdb_idx}"
        pdb_labels.append(pdb_label)
        row = rows_by_idx.get(idx)
        ch_label = row.get("chothia") if row else ""
        chothia_labels.append(ch_label or pdb_label)
    return pdb_labels, chothia_labels


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
    pdb_a: Path,
    chain_a: str,
    pdb_b: Path,
    chains_b: List[str],
    numbering_path: Path | None,
    cutoff: float,
    buffer: int,
    output: Path,
    title: str | None,
):
    parser = PDBParser(QUIET=True)
    structure_a = parser.get_structure("A", pdb_a)
    structure_b = parser.get_structure("B", pdb_b)

    chain_a_obj = structure_a[0][chain_a]
    coords_a, labels_a_raw = _iter_ca(chain_a_obj)
    labels_a = [f"{chain_a}{pdb_idx}{icode}" if icode else f"{chain_a}{pdb_idx}" for pdb_idx, icode, _ in labels_a_raw]

    num_map = _load_numbering_map(numbering_path) if numbering_path else {}

    coords_b_list: List[np.ndarray] = []
    labels_b_pdb: List[str] = []
    labels_b_chothia: List[str] = []
    for cid in chains_b:
        chain = structure_b[0][cid]
        coords_b, labels_raw_b = _iter_ca(chain)
        pdb_labels, ch_labels = _labels_for_chain(cid, labels_raw_b, num_map)
        coords_b_list.append(coords_b)
        labels_b_pdb.extend(pdb_labels)
        labels_b_chothia.extend(ch_labels)
    coords_b_concat = np.vstack(coords_b_list)

    distances = _compute_contacts(coords_a, coords_b_concat)
    cropped, (row_idx, col_idx) = _crop_matrix(distances, cutoff, buffer)
    labels_a_crop = [labels_a[i] for i in row_idx]
    labels_b_pdb_crop = [labels_b_pdb[j] for j in col_idx]
    labels_b_chothia_crop = [labels_b_chothia[j] for j in col_idx]

    masked = np.where(cropped <= cutoff, cropped, np.nan)
    cmap = plt.cm.coolwarm
    cmap.set_bad("white")

    fig_w = max(6.0, len(labels_b_pdb_crop) * 0.25)
    fig_h = max(5.0, len(labels_a_crop) * 0.25)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(masked, origin="lower", cmap=cmap, aspect="auto")
    cbar = fig.colorbar(im, ax=ax, label="CA-CA distance (Å)")
    cbar.ax.tick_params(labelsize=9)

    ax.set_xticks(np.arange(len(labels_b_pdb_crop)))
    ax.set_yticks(np.arange(len(labels_a_crop)))
    ax.set_xticklabels(labels_b_pdb_crop, rotation=90, fontsize=8)
    ax.set_yticklabels(labels_a_crop, fontsize=8)
    ax.set_xlabel("Residue (PDB index, chain B)")
    ax.set_ylabel("Residue (PDB index, chain A)")

    ax_top = ax.secondary_xaxis('top')
    ax_top.set_xticks(np.arange(len(labels_b_chothia_crop)))
    ax_top.set_xticklabels(labels_b_chothia_crop, rotation=90, fontsize=8)
    ax_top.set_xlabel("Chothia numbering (chain B)")

    name_a = pdb_a.name
    name_b = pdb_b.name
    ax.set_title(title or f"Contact map: {name_a} vs {name_b}")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    print(f"Saved contact map to {output}")


def main(argv: Sequence[str] | None = None):
    ap = argparse.ArgumentParser(description="Plot CA contact map with dual axes (PDB + Chothia)")
    ap.add_argument("--pdb-a", required=True, type=Path, help="PDB for chain A")
    ap.add_argument("--chain-a", required=True, help="Chain ID for protein A")
    ap.add_argument("--pdb-b", required=True, type=Path, help="PDB for protein/antibody B")
    ap.add_argument("--chains-b", required=True, help="Chain IDs for protein/antibody B (space/comma)")
    ap.add_argument("--numbering-map", type=Path, default=None, help="CSV/JSON numbering map (from build_antibody_numbering.py)")
    ap.add_argument("--cutoff", type=float, default=8.0, help="CA-CA distance cutoff (Å)")
    ap.add_argument("--buffer", type=int, default=2, help="Residue buffer around contacting region")
    ap.add_argument("--output", type=Path, default=Path("contact_map.png"), help="Output image path")
    ap.add_argument("--title", type=str, default=None, help="Plot title override")
    args = ap.parse_args(argv)

    chains_b = _clean_list(args.chains_b)
    if not chains_b:
        raise SystemExit("No chains provided for B")

    plot_contact_map(
        pdb_a=args.pdb_a,
        chain_a=args.chain_a,
        pdb_b=args.pdb_b,
        chains_b=chains_b,
        numbering_path=args.numbering_map,
        cutoff=args.cutoff,
        buffer=args.buffer,
        output=args.output,
        title=args.title,
    )


if __name__ == "__main__":
    main()
