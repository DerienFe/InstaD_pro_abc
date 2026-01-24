#!/usr/bin/env python3
"""Generate PDB↔Chothia numbering maps (CSV/JSON) with optional CDR annotation.

- Uses ANARCI to renumber provided antibody chains (Chothia scheme by default).
- Emits per-residue mapping: chain, global_idx, chain_idx, pdb_idx, icode, chothia_label, aa, cdr_label.
- JSON includes CDR summaries with sequences and indices.
"""
import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import PDB
from Bio.SeqUtils import seq1

# Chothia CDR boundaries (inclusive) in numbering space
CDR_BOUNDARIES = {
    "H": {
        "H1": (26, 32),
        "H2": (52, 56),
        "H3": (95, 102),
    },
    "L": {
        "L1": (24, 34),
        "L2": (50, 56),
        "L3": (89, 97),
    },
}


def _clean_list(raw: str) -> List[str]:
    if not raw:
        return []
    parts = raw.replace(",", " ").split()
    uniq: List[str] = []
    for p in parts:
        if p not in uniq:
            uniq.append(p)
    return uniq


def _iter_residues(chain) -> Tuple[List, List[str]]:
    residues = []
    aas: List[str] = []
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
        aas.append(aa)
    return residues, aas


def _run_anarci(seq: str, scheme: str) -> List[Optional[str]]:
    from anarci import anarci  # type: ignore

    # Ensure hmmscan from the active env is discoverable
    hmmscan = shutil.which("hmmscan")
    if not hmmscan:
        candidate = Path(sys.executable).parent / "hmmscan"
        if candidate.exists():
            os.environ["PATH"] = f"{candidate.parent}:{os.environ.get('PATH', '')}"

    numbering, _, _ = anarci([("chain", seq)], scheme=scheme)
    if not numbering or numbering[0] is None or numbering[0][0] is None:
        raise ValueError("ANARCI could not number the provided sequence")

    domain_numbering, start_idx, _ = numbering[0][0]
    labels: List[Optional[str]] = []
    for entry in domain_numbering:
        if entry is None:
            labels.append(None)
            continue
        (pos, ins), _aa = entry
        label = f"{pos}{ins.strip()}" if ins and ins.strip() else str(pos)
        labels.append(label)
    labels = labels[start_idx:]
    if len(labels) < len(seq):
        labels.extend([None] * (len(seq) - len(labels)))
    if len(labels) > len(seq):
        labels = labels[:len(seq)]
    return labels


def _cdr_label(chain_kind: Optional[str], label: Optional[str]) -> str:
    if not chain_kind or not label or chain_kind not in CDR_BOUNDARIES:
        return ""
    try:
        base_num = int("".join(ch for ch in label if ch.isdigit()))
    except ValueError:
        return ""
    for name, (lo, hi) in CDR_BOUNDARIES[chain_kind].items():
        if lo <= base_num <= hi:
            return name
    return ""


def build_maps(pdb_path: Path, chain_ids: List[str], heavy_ids: List[str], light_ids: List[str], scheme: str) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, pdb_path)
    available = {c.id for c in structure.get_chains()}
    missing = [c for c in chain_ids if c not in available]
    if missing:
        raise ValueError(f"Chains not found in {pdb_path.name}: {missing}")

    rows: List[Dict[str, object]] = []
    chains_json: Dict[str, object] = {}
    global_idx = 1
    for cid in chain_ids:
        chain = structure[0][cid]
        residues, aas = _iter_residues(chain)
        if not residues:
            raise ValueError(f"No standard residues with CA atoms found in chain {cid}")
        chain_kind: Optional[str]
        if cid in heavy_ids:
            chain_kind = "H"
        elif cid in light_ids:
            chain_kind = "L"
        else:
            chain_kind = None

        chothia_labels = _run_anarci("".join(aas), scheme)
        if len(chothia_labels) != len(residues):
            print(
                f"[warn] Chain {cid}: ANARCI labels {len(chothia_labels)} vs residues {len(residues)}; padding/truncating",
                file=sys.stderr,
            )
            if len(chothia_labels) < len(residues):
                chothia_labels.extend([None] * (len(residues) - len(chothia_labels)))
            else:
                chothia_labels = chothia_labels[: len(residues)]

        chain_rows: List[Dict[str, object]] = []
        cdr_bins: Dict[str, Dict[str, object]] = {}
        for idx, (res, aa, ch_label) in enumerate(zip(residues, aas, chothia_labels), start=1):
            pdb_idx = int(res.id[1])
            icode = res.id[2].strip() or ""
            cdr = _cdr_label(chain_kind, ch_label)
            entry = {
                "chain": cid,
                "global_idx": global_idx,
                "chain_idx": idx,
                "pdb_idx": pdb_idx,
                "icode": icode,
                "chothia": ch_label or "",
                "aa": aa,
                "cdr": cdr,
            }
            rows.append(entry)
            chain_rows.append(entry)
            global_idx += 1
            if cdr:
                bucket = cdr_bins.setdefault(cdr, {"chain_idx": [], "pdb_idx": [], "chothia": [], "seq": []})
                bucket["chain_idx"].append(idx)
                bucket["pdb_idx"].append(pdb_idx)
                bucket["chothia"].append(ch_label or "")
                bucket["seq"].append(aa)

        # finalize CDR seq strings
        for cdr_name, bucket in cdr_bins.items():
            bucket["sequence"] = "".join(bucket.pop("seq"))

        chains_json[cid] = {
            "length": len(chain_rows),
            "rows": chain_rows,
            "cdrs": cdr_bins,
            "type": chain_kind or "unknown",
        }

    json_out = {
        "pdb": str(pdb_path),
        "scheme": scheme,
        "chains": chain_ids,
        "map": chains_json,
    }
    return rows, json_out


def write_csv(rows: List[Dict[str, object]], path: Path) -> None:
    header = ["chain", "global_idx", "chain_idx", "pdb_idx", "icode", "chothia", "aa", "cdr"]
    lines = [",".join(header)]
    for r in rows:
        lines.append(
            f"{r['chain']},{r['global_idx']},{r['chain_idx']},{r['pdb_idx']},{r['icode']},{r['chothia']},{r['aa']},{r['cdr']}"
        )
    path.write_text("\n".join(lines) + "\n")


def write_json(payload: Dict[str, object], path: Path) -> None:
    path.write_text(json.dumps(payload, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(description="Build Chothia numbering maps for antibody chains")
    ap.add_argument("--pdb", required=True, type=Path, help="PDB file containing the antibody chains")
    ap.add_argument("--chains", required=True, help="Chain IDs to number (space or comma separated)")
    ap.add_argument("--heavy", default="H", help="Heavy chain IDs (space/comma). Default: H")
    ap.add_argument("--light", default="L", help="Light chain IDs (space/comma). Default: L")
    ap.add_argument("--scheme", default="chothia", help="Numbering scheme for ANARCI (default: chothia)")
    ap.add_argument("--out-csv", type=Path, default=None, help="Output CSV path")
    ap.add_argument("--out-json", type=Path, default=None, help="Output JSON path")
    args = ap.parse_args()

    chains = _clean_list(args.chains)
    if not chains:
        raise SystemExit("No chains provided")
    heavy_ids = _clean_list(args.heavy)
    light_ids = _clean_list(args.light)

    rows, payload = build_maps(args.pdb, chains, heavy_ids, light_ids, args.scheme)

    stem = args.pdb.stem
    out_csv = args.out_csv or args.pdb.with_name(f"{stem}_chothia_map.csv")
    out_json = args.out_json or args.pdb.with_name(f"{stem}_chothia_map.json")
    write_csv(rows, out_csv)
    write_json(payload, out_json)
    print(f"Wrote numbering map CSV to {out_csv}")
    print(f"Wrote numbering map JSON to {out_json}")


if __name__ == "__main__":
    main()
