#!/usr/bin/env python3
"""Compute per-residue NLL deltas between bound/unbound runs and plot.

- Inputs: bound/unbound design_nll CSVs (columns: chain,pdb_idx,chain_idx,global_idx,nll)
- Output: CSV with deltas, PNG plot, and top-5 increase/decrease text summary.
- Optional: Chothia map (JSON from build_antibody_numbering.py) to annotate ticks for antibody chains.
"""
import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


def load_design_nll(path: Path) -> Dict[Tuple[str, int], Dict[str, float | str | None]]:
    table: Dict[Tuple[str, int], Dict[str, float | str | None]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            chain = row.get('chain') or row.get('Chain')
            if chain is None:
                continue
            pdb_idx = row.get('pdb_idx') or row.get('pdb') or row.get('pdb_index')
            if pdb_idx is None or row['nll'] is None:
                continue
            try:
                pdb_idx_int = int(pdb_idx)
                nll_val = float(row['nll'])
            except ValueError:
                continue
            chain_idx = int(row['chain_idx']) if row.get('chain_idx') else None
            aa_val = row.get('aa') if row.get('aa') else None
            table[(chain, pdb_idx_int)] = {
                'nll': nll_val,
                'chain_idx': chain_idx,
                'global_idx': int(row['global_idx']) if row.get('global_idx') else None,
                'aa': aa_val,
            }
    return table


def load_chothia_map(path: Path | None) -> Dict[Tuple[str, int], Dict[str, str]]:
    """Return per-residue annotations: chothia label and one-letter AA."""
    if path is None or not path.exists():
        return {}
    data = json.loads(path.read_text())
    out: Dict[Tuple[str, int], Dict[str, str]] = {}
    chains = data.get('map', {}) if isinstance(data, dict) else {}
    for cid, payload in chains.items():
        rows = payload.get('rows', []) if isinstance(payload, dict) else []
        for r in rows:
            try:
                pdb_idx = int(r.get('pdb_idx', 0))
            except (TypeError, ValueError):
                continue
            label = r.get('chothia', '') or ''
            aa = (r.get('aa') or '').strip()
            out[(cid, pdb_idx)] = {'chothia': label, 'aa': aa}
    return out


def load_generic_map(path: Path | None) -> Dict[Tuple[str, int], str]:
    """Load generic pdb-index -> AA map produced by preprocess_chains (_map.json)."""
    if path is None or not path.exists():
        return {}
    data = json.loads(path.read_text())
    out: Dict[Tuple[str, int], str] = {}
    chains = data.get('map', {}) if isinstance(data, dict) else {}
    for cid, payload in chains.items():
        rows = payload.get('rows', []) if isinstance(payload, dict) else []
        for r in rows:
            try:
                pdb_idx = int(r.get('pdb_idx', 0))
            except (TypeError, ValueError):
                continue
            aa = (r.get('aa') or '').strip()
            if aa:
                out[(cid, pdb_idx)] = aa
    return out


def compute_delta(bound: Dict[Tuple[str, int], Dict[str, float | str | None]], unbound: Dict[Tuple[str, int], Dict[str, float | str | None]]):
    deltas = []
    keys = set(bound.keys()) & set(unbound.keys())
    for key in sorted(keys):
        b = bound[key]['nll']
        u = unbound[key]['nll']
        deltas.append({
            'chain': key[0],
            'pdb_idx': key[1],
            'chain_idx': bound[key].get('chain_idx'),
            'global_idx': bound[key].get('global_idx'),
            'bound_nll': b,
            'unbound_nll': u,
            'delta': b - u,
            'aa': bound[key].get('aa') or unbound[key].get('aa'),
        })
    return deltas


def write_delta_csv(deltas: List[Dict[str, object]], path: Path, chothia: Dict[Tuple[str, int], Dict[str, str]], pdb_map: Dict[Tuple[str, int], str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ['chain', 'pdb_idx', 'chain_idx', 'global_idx', 'bound_nll', 'unbound_nll', 'delta', 'aa', 'chothia']
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        for d in deltas:
            key = (d['chain'], d['pdb_idx'])
            aa_val = d.get('aa') or chothia.get(key, {}).get('aa', '') or pdb_map.get(key, '')
            writer.writerow({
                **d,
                'aa': aa_val,
                'chothia': chothia.get(key, {}).get('chothia', ''),
            })


def top_changes(deltas: List[Dict[str, object]], chothia: Dict[Tuple[str, int], Dict[str, str]], pdb_map: Dict[Tuple[str, int], str], k: int = 5):
    inc = sorted(deltas, key=lambda x: x['delta'], reverse=True)[:k]
    dec = sorted(deltas, key=lambda x: x['delta'])[:k]
    def fmt(lst):
        lines = []
        for d in lst:
            key = (d['chain'], d['pdb_idx'])
            ann = chothia.get(key, {})
            ch = ann.get('chothia', '')
            aa = ann.get('aa', '') or pdb_map.get(key, '') or d.get('aa', '')
            aa_tag = f"{aa}{d['pdb_idx']}" if aa else str(d['pdb_idx'])
            ch_tag = f" (chothia {ch})" if ch else ''
            lines.append(f"{d['chain']} {aa_tag}{ch_tag}: delta={d['delta']:.3f} bound={d['bound_nll']:.3f} unbound={d['unbound_nll']:.3f}")
        return lines
    return fmt(inc), fmt(dec)


def plot_deltas(deltas: List[Dict[str, object]], chothia: Dict[Tuple[str, int], Dict[str, str]], pdb_map: Dict[Tuple[str, int], str], title: str, out_png: Path):
    chains = sorted({d['chain'] for d in deltas})
    fig, axes = plt.subplots(len(chains), 1, figsize=(10, max(3, 3 * len(chains))), sharey=True)
    if hasattr(axes, 'flatten'):
        axes = list(axes.flatten())
    elif not isinstance(axes, (list, tuple)):
        axes = [axes]
    for ax, chain in zip(axes, chains):
        subset = sorted([d for d in deltas if d['chain'] == chain], key=lambda x: x['pdb_idx'])
        raw_idx = [d['pdb_idx'] for d in subset]
        y = [d['delta'] for d in subset]

        # Compress gaps between disjoint regions to a small spacer to avoid large empty spans on the x-axis.
        gap_cap = 5  # maximum virtual spacer length for large jumps
        x = []
        cursor = 0
        prev = None
        for idx in raw_idx:
            if prev is None:
                cursor = 0
            else:
                gap = idx - prev
                spacer = 1 if gap <= 1 else min(gap - 1, gap_cap) + 1
                cursor += spacer
            x.append(cursor)
            prev = idx

        ax.axhline(0.0, color='gray', lw=1, linestyle='--')
        ax.bar(x, y, width=0.8, align='center', color='#4C72B0')
        # Force one tick per residue in the designable set.
        ax.set_xticks(x)
        labels = []
        for pdb_idx, d in zip(raw_idx, subset):
            ann = chothia.get((chain, pdb_idx), {})
            aa = ann.get('aa', '') or pdb_map.get((chain, pdb_idx), '') or d.get('aa', '') or ''
            labels.append(f"{aa}{pdb_idx}" if aa else str(pdb_idx))
        ax.set_xticklabels(labels)
        ax.set_xlim(min(x) - 0.6, max(x) + 0.6)
        ax.set_xlabel(f"Chain {chain} PDB index (compressed gaps)")
        ax.set_ylabel('Delta NLL (bound - unbound)')
        ax.set_title(f"{title} chain {chain}")
        # Top axis with Chothia labels if available
        ticks = x
        top_labels = [chothia.get((chain, pdb_idx), {}).get('chothia', '') for pdb_idx in raw_idx]
        if any(top_labels):
            ax_top = ax.secondary_xaxis('top')
            ax_top.set_xticks(ticks)
            ax_top.set_xticklabels(top_labels, rotation=90, fontsize=8)
            ax_top.set_xlabel('Chothia index')
        ax.tick_params(axis='x', rotation=90, labelsize=8)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)


def main():
    ap = argparse.ArgumentParser(description="Analyze bound vs unbound design_nll deltas and plot")
    ap.add_argument('--bound-csv', required=True, type=Path, help='design_nll CSV for bound state')
    ap.add_argument('--unbound-csv', required=True, type=Path, help='design_nll CSV for unbound state')
    ap.add_argument('--label', required=True, help='Label for outputs (e.g., 7z0x or 6m0j)')
    ap.add_argument('--chothia-map', type=Path, default=None, help='Optional Chothia map JSON for antibody chains')
    ap.add_argument('--map-json', type=Path, default=None, help='Optional generic map JSON from preprocess_chains (_map.json)')
    ap.add_argument('--out-dir', type=Path, default=Path('proteinmpnn_run/results/analyzed'))
    args = ap.parse_args()

    bound = load_design_nll(args.bound_csv)
    unbound = load_design_nll(args.unbound_csv)
    chothia = load_chothia_map(args.chothia_map)
    pdb_map = load_generic_map(args.map_json)

    deltas = compute_delta(bound, unbound)
    if not deltas:
        raise SystemExit('No overlapping residues between bound and unbound inputs')

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    delta_csv = out_dir / f"{args.label}_delta.csv"
    delta_png = out_dir / f"{args.label}_delta.png"
    delta_txt = out_dir / f"{args.label}_top.txt"

    write_delta_csv(deltas, delta_csv, chothia, pdb_map)
    plot_deltas(deltas, chothia, pdb_map, args.label, delta_png)

    inc, dec = top_changes(deltas, chothia, pdb_map, k=5)
    delta_txt.write_text(
        "Top +5 (bound worse):\n" + "\n".join(inc) + "\n\n" + "Top -5 (bound better):\n" + "\n".join(dec) + "\n"
    )
    print(f"Wrote delta CSV to {delta_csv}")
    print(f"Wrote plot to {delta_png}")
    print(f"Wrote top changes to {delta_txt}")


if __name__ == '__main__':
    main()
