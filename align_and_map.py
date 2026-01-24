#!/usr/bin/env python3
import argparse
import json
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from anarci import run_anarci as anarci_run
except ImportError:  # pragma: no cover - optional dependency
    anarci_run = None

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z"
}

CDR_RANGES = {
    "H": {
        "CDR-H1": (26, 32),
        "CDR-H2": (52, 56),
        "CDR-H3": (95, 102),
    },
    "L": {
        "CDR-L1": (24, 34),
        "CDR-L2": (50, 56),
        "CDR-L3": (89, 97),
    },
}


@dataclass
class ResidueRecord:
    aa: str
    resname: str
    auth_seq_id: Optional[str]
    label_seq_id: Optional[str]
    ins_code: Optional[str]
    chain_id: str


def three_to_one(resname: str) -> str:
    return THREE_TO_ONE.get(resname.upper(), "X")


def read_fasta(path: Path) -> Tuple[str, str]:
    header = None
    seq_lines: List[str] = []
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].strip() or path.stem
            else:
                seq_lines.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA header in {path}")
    return header, "".join(seq_lines).replace(" ", "").upper()


def needleman_wunsch(seq_a: str, seq_b: str) -> Tuple[str, str]:
    match_score = 2
    mismatch_score = -1
    gap_penalty = -2
    len_a, len_b = len(seq_a), len(seq_b)
    score = [[0] * (len_b + 1) for _ in range(len_a + 1)]
    pointer = [[None] * (len_b + 1) for _ in range(len_a + 1)]

    for i in range(1, len_a + 1):
        score[i][0] = i * gap_penalty
        pointer[i][0] = "up"
    for j in range(1, len_b + 1):
        score[0][j] = j * gap_penalty
        pointer[0][j] = "left"

    for i in range(1, len_a + 1):
        for j in range(1, len_b + 1):
            diag = score[i - 1][j - 1] + (match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_score)
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                pointer[i][j] = "diag"
            elif best == up:
                pointer[i][j] = "up"
            else:
                pointer[i][j] = "left"

    align_a: List[str] = []
    align_b: List[str] = []
    i, j = len_a, len_b
    while i > 0 or j > 0:
        move = pointer[i][j]
        if move == "diag":
            align_a.append(seq_a[i - 1])
            align_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            align_a.append(seq_a[i - 1])
            align_b.append("-")
            i -= 1
        else:
            align_a.append("-")
            align_b.append(seq_b[j - 1])
            j -= 1
    return "".join(reversed(align_a)), "".join(reversed(align_b))


def _mmcif_value(fields: List[str], indices: Dict[str, int], names: Sequence[str]) -> str:
    for name in names:
        idx = indices.get(name)
        if idx is not None and idx < len(fields):
            value = fields[idx]
            if value not in {"?", "."}:
                return value
    return ""


def parse_mmcif_chain(structure_path: Path, chain_id: str) -> List[ResidueRecord]:
    headers: List[str] = []
    indices: Dict[str, int] = {}
    within_atom_loop = False
    last_key: Tuple[Optional[str], Optional[str], Optional[str]] = (None, None, None)
    target_chain = chain_id.strip()
    residues: List[ResidueRecord] = []

    with structure_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                within_atom_loop = False
                headers = []
                indices = {}
                continue
            if stripped.startswith("loop_"):
                continue
            if stripped.startswith("_"):
                if stripped.startswith("_atom_site."):
                    headers.append(stripped)
                continue
            if not headers:
                continue
            if not within_atom_loop:
                within_atom_loop = True
                indices = {name: idx for idx, name in enumerate(headers)}
            fields = shlex.split(stripped)
            if len(fields) < len(headers):
                continue
            group = _mmcif_value(fields, indices, ["_atom_site.group_PDB"]).upper()
            if group not in {"ATOM", "HETATM"}:
                continue
            auth_chain = _mmcif_value(
                fields,
                indices,
                ["_atom_site.auth_asym_id", "_atom_site.label_asym_id"],
            ).strip()
            if auth_chain != target_chain:
                continue
            resname = _mmcif_value(
                fields,
                indices,
                ["_atom_site.auth_comp_id", "_atom_site.label_comp_id"],
            ).strip()
            auth_seq_id = _mmcif_value(
                fields,
                indices,
                ["_atom_site.auth_seq_id"],
            ).strip() or None
            label_seq_id = _mmcif_value(
                fields,
                indices,
                ["_atom_site.label_seq_id"],
            ).strip() or None
            ins_code = _mmcif_value(fields, indices, ["_atom_site.pdbx_PDB_ins_code"]).strip() or None
            current_key = (auth_chain, auth_seq_id, ins_code)
            if current_key == last_key:
                continue
            last_key = current_key
            residues.append(
                ResidueRecord(
                    aa=three_to_one(resname),
                    resname=resname,
                    auth_seq_id=auth_seq_id,
                    label_seq_id=label_seq_id,
                    ins_code=ins_code,
                    chain_id=auth_chain,
                )
            )
    if not residues:
        raise ValueError(f"No residues parsed for chain {chain_id} in {structure_path}")
    return residues


def build_alignment_mapping(
    align_ref: str,
    align_query: str,
    ref_records: List[ResidueRecord],
    query_records: List[ResidueRecord],
) -> Tuple[List[Dict[str, Optional[str]]], str, str]:
    if len(align_ref) != len(align_query):
        raise ValueError("Aligned sequences must have equal length")
    ref_idx = -1
    query_idx = -1
    mapping: List[Dict[str, Optional[str]]] = []
    coverage_ref: List[str] = []
    coverage_query: List[str] = []

    for ref_char, query_char in zip(align_ref, align_query):
        ref_meta: Optional[ResidueRecord] = None
        query_meta: Optional[ResidueRecord] = None
        if ref_char != "-":
            ref_idx += 1
            if ref_idx >= len(ref_records):
                raise IndexError("Reference index exceeds residue annotations")
            ref_meta = ref_records[ref_idx]
        if query_char != "-":
            query_idx += 1
            if query_idx >= len(query_records):
                raise IndexError("Query index exceeds residue annotations")
            query_meta = query_records[query_idx]
        if ref_char != "-" and query_char != "-":
            coverage_ref.append(ref_char)
            coverage_query.append(query_char)
            mapping.append(
                {
                    "reference_position": ref_idx,
                    "reference_residue": ref_meta.aa if ref_meta else None,
                    "reference_auth_seq_id": ref_meta.auth_seq_id if ref_meta else None,
                    "reference_label_seq_id": ref_meta.label_seq_id if ref_meta else None,
                    "reference_ins_code": ref_meta.ins_code if ref_meta else None,
                    "reference_chain": ref_meta.chain_id if ref_meta else None,
                    "query_position": query_idx,
                    "query_residue": query_meta.aa if query_meta else None,
                    "query_auth_seq_id": query_meta.auth_seq_id if query_meta else None,
                    "query_label_seq_id": query_meta.label_seq_id if query_meta else None,
                    "query_ins_code": query_meta.ins_code if query_meta else None,
                    "query_chain": query_meta.chain_id if query_meta else None,
                }
            )
    return mapping, "".join(coverage_ref), "".join(coverage_query)


def write_fasta(output_path: Path, entries: List[Tuple[str, str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="ascii") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def write_alignment_json(
    output_path: Path,
    reference_info: Dict[str, str],
    query_info: Dict[str, str],
    mapping: List[Dict[str, Optional[str]]],
) -> None:
    payload = {
        "reference": reference_info,
        "query": query_info,
        "aligned_pairs": mapping,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def run_anarci_mapping(
    fasta_path: Path,
    chain_label: str,
    chain_type_hint: str,
    output_path: Path,
) -> None:
    if anarci_run is None:
        print("[warn] ANARCI is not installed; skipping antibody numbering", file=sys.stderr)
        return
    if not fasta_path.exists():
        print(f"[warn] Missing ANARCI FASTA: {fasta_path}", file=sys.stderr)
        return
    header, sequence = read_fasta(fasta_path)
    result = anarci_run([(header, sequence)], scheme="chothia")
    numbering_payload = _extract_anarci_numbering(result, chain_type_hint)
    if numbering_payload is None:
        print("[warn] Unable to extract ANARCI numbering", file=sys.stderr)
        return
    chain_type, numbering = numbering_payload
    numbered_seq = "".join(residue for _, residue in numbering)
    cdr_ranges = CDR_RANGES.get(chain_type, {})
    cdr_sequences: Dict[str, List[str]] = {name: [] for name in cdr_ranges}
    bin_mask: List[str] = []
    for (position, insertion), residue in numbering:
        base_pos = position if isinstance(position, int) else int(position)
        in_cdr = False
        for cdr_name, (start, end) in cdr_ranges.items():
            if start <= base_pos <= end:
                in_cdr = True
                cdr_sequences[cdr_name].append(residue)
        bin_mask.append("1" if in_cdr else "0")
    mapping_entry: Dict[str, str] = {
        "chain_label": chain_label,
        "chain_type": chain_type,
        "original_seq": sequence,
        "numbered_seq": numbered_seq,
        "cdr_bin_mask": "".join(bin_mask),
    }
    for cdr_name, chars in cdr_sequences.items():
        mapping_entry[cdr_name] = "".join(chars)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump([mapping_entry], handle, indent=2)
    print(f"[info] Wrote ANARCI mapping to {output_path}")


def _extract_anarci_numbering(result, chain_type_hint: str) -> Optional[Tuple[str, List[Tuple[Tuple[int, str], str]]]]:
    domain_results = result[0]
    if not domain_results:
        return None
    first_entry = domain_results[0]
    if isinstance(first_entry, tuple) and len(first_entry) >= 2:
        domains = first_entry[1]
    else:
        domains = first_entry
    if not domains:
        return None
    chosen = None
    for domain in domains:
        chain_type = None
        numbering = None
        if isinstance(domain, dict):
            chain_type = domain.get("chain_type")
            numbering = domain.get("numbering")
        elif isinstance(domain, tuple):
            if len(domain) >= 3 and isinstance(domain[2], list):
                chain_type = domain[0]
                numbering = domain[2]
            elif len(domain) >= 4 and isinstance(domain[3], list):
                chain_type = domain[0]
                numbering = domain[3]
        if numbering is None:
            continue
        if chain_type_hint and chain_type and chain_type.upper() != chain_type_hint.upper():
            continue
        chosen = (chain_type or chain_type_hint, numbering)
        break
    if chosen is None:
        first_domain = domains[0]
        if isinstance(first_domain, dict):
            return first_domain.get("chain_type", chain_type_hint), first_domain.get("numbering", [])
        if isinstance(first_domain, tuple):
            chain_type = first_domain[0] if first_domain else chain_type_hint
            numbering = first_domain[2] if len(first_domain) > 2 else []
            return chain_type, numbering
        return None
    formatted_numbering: List[Tuple[Tuple[int, str], str]] = []
    chain_type, numbering = chosen
    for item in numbering:
        if not item:
            continue
        if isinstance(item[0], tuple):
            position_tuple = item[0]
            residue = item[1] if len(item) > 1 else "X"
        else:
            position_tuple = item[:2] if len(item) >= 2 else (item[0], "")
            residue = item[2] if len(item) >= 3 else "X"
        pos_index = position_tuple[0]
        insert_code = position_tuple[1] if len(position_tuple) > 1 else ""
        formatted_numbering.append(((int(pos_index), str(insert_code)), residue))
    return chain_type or chain_type_hint, formatted_numbering


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Align FASTA sequences and build residue mappings.")
    parser.add_argument("--reference-fasta", type=Path, default=Path("data/6M0J_E.fasta"))
    parser.add_argument("--reference-structure", type=Path, default=Path("data/6M0J.cif"))
    parser.add_argument("--reference-chain", default="E")
    parser.add_argument("--query-fasta", type=Path, default=Path("data/7Z0X_R.fasta"))
    parser.add_argument("--query-structure", type=Path, default=Path("data/7Z0X.cif"))
    parser.add_argument("--query-chain", default="R")
    parser.add_argument("--output-dir", type=Path, default=Path("data"))
    parser.add_argument("--coverage-fasta-name", default="max_coverage_alignment.fasta")
    parser.add_argument("--mapping-json-name", default="residue_alignment_mapping.json")
    parser.add_argument("--anarci-fasta", type=Path, default=Path("data/7Z0X_H.fasta"))
    parser.add_argument("--anarci-chain-label", default="7Z0X_H")
    parser.add_argument("--anarci-chain-type", choices=["H", "L"], default="H")
    parser.add_argument("--anarci-json-name", default="anarci_chothia_mapping.json")
    parser.add_argument("--skip-anarci", action="store_true")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    for path in [args.reference_fasta, args.query_fasta, args.reference_structure, args.query_structure]:
        if not path.exists():
            parser.error(f"Required file not found: {path}")

    ref_header, ref_seq = read_fasta(args.reference_fasta)
    query_header, query_seq = read_fasta(args.query_fasta)
    ref_records = parse_mmcif_chain(args.reference_structure, args.reference_chain)
    query_records = parse_mmcif_chain(args.query_structure, args.query_chain)

    ref_struct_seq = "".join(res.aa for res in ref_records)
    query_struct_seq = "".join(res.aa for res in query_records)
    if len(ref_struct_seq) != len(ref_records):
        parser.error("Reference structure parsing failed to map residues")
    if len(query_struct_seq) != len(query_records):
        parser.error("Query structure parsing failed to map residues")

    align_ref, align_query = needleman_wunsch(ref_seq, query_seq)
    mapping, coverage_ref, coverage_query = build_alignment_mapping(align_ref, align_query, ref_records, query_records)

    output_dir = args.output_dir
    coverage_path = output_dir / args.coverage_fasta_name
    mapping_path = output_dir / args.mapping_json_name

    write_fasta(
        coverage_path,
        [
            (f"{ref_header}_maxcov", coverage_ref),
            (f"{query_header}_maxcov", coverage_query),
        ],
    )
    print(f"[info] Wrote coverage FASTA to {coverage_path}")

    reference_info = {
        "label": ref_header,
        "chain": args.reference_chain,
        "structure": str(args.reference_structure),
        "sequence_length": len(ref_seq),
    }
    query_info = {
        "label": query_header,
        "chain": args.query_chain,
        "structure": str(args.query_structure),
        "sequence_length": len(query_seq),
    }
    write_alignment_json(mapping_path, reference_info, query_info, mapping)
    print(f"[info] Wrote residue mapping JSON to {mapping_path}")

    if not args.skip_anarci:
        anarci_output_path = output_dir / args.anarci_json_name
        run_anarci_mapping(args.anarci_fasta, args.anarci_chain_label, args.anarci_chain_type, anarci_output_path)


if __name__ == "__main__":
    main()
