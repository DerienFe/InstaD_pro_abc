import argparse
import shlex
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z"
}


def three_to_one(resname: str) -> str:
    """Return the one-letter code for a three-letter residue name."""
    return THREE_TO_ONE.get(resname.upper(), "X")


def parse_args() -> Tuple[argparse.Namespace, argparse.ArgumentParser]:
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences from PDB or mmCIF files."
    )
    parser.add_argument(
        "--structure",
        type=Path,
        help="Path to the structure file (PDB or mmCIF). Defaults to legacy pro/ batch when omitted.",
    )
    parser.add_argument(
        "-c",
        "--chain",
        dest="chains",
        action="append",
        help="Chain identifier to extract. Repeat for multiple chains.",
    )
    parser.add_argument(
        "--label",
        help="Identifier used for FASTA headers and filenames. Defaults to the structure stem.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where FASTA files are written (defaults to CWD).",
    )
    return parser.parse_args(), parser


def parse_pdb_chains(pdb_path: Path) -> Dict[str, List[str]]:
    chains: Dict[str, List[str]] = {}
    last_key = (None, None, None)

    with pdb_path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            resname = line[17:20].strip()
            chain_id = line[21].strip() or "_"
            res_seq = line[22:26].strip()
            icode = line[26].strip()

            current_key = (chain_id, res_seq, icode)
            if current_key == last_key:
                continue
            last_key = current_key

            chains.setdefault(chain_id, []).append(three_to_one(resname))

    return chains


def _mmcif_value(fields: List[str], indices: Dict[str, int], names: Sequence[str]) -> str:
    for name in names:
        idx = indices.get(name)
        if idx is not None and idx < len(fields):
            value = fields[idx]
            if value not in {".", "?"}:
                return value
    return ""


def parse_mmcif_chains(cif_path: Path) -> Dict[str, List[str]]:
    chains: Dict[str, List[str]] = {}
    headers: List[str] = []
    indices: Dict[str, int] = {}
    within_atom_loop = False
    last_key = (None, None, None)

    with cif_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if within_atom_loop:
                    within_atom_loop = False
                    headers = []
                    indices = {}
                continue
            if stripped.startswith("loop_"):
                if within_atom_loop:
                    within_atom_loop = False
                    headers = []
                    indices = {}
                continue
            if stripped.startswith("_"):
                if stripped.startswith("_atom_site."):
                    headers.append(stripped)
                else:
                    if within_atom_loop:
                        within_atom_loop = False
                        headers = []
                        indices = {}
                continue
            if not headers:
                continue
            if not within_atom_loop:
                within_atom_loop = True
                indices = {name: idx for idx, name in enumerate(headers)}
            fields = shlex.split(stripped)
            if len(fields) < len(headers):
                continue
            group = _mmcif_value(fields, indices, ["_atom_site.group_PDB"])
            if group and group.upper() not in {"ATOM", "HETATM"}:
                continue
            resname = _mmcif_value(
                fields,
                indices,
                ["_atom_site.auth_comp_id", "_atom_site.label_comp_id"],
            )
            chain_id = _mmcif_value(
                fields,
                indices,
                ["_atom_site.auth_asym_id", "_atom_site.label_asym_id"],
            ).strip() or "_"
            res_seq = _mmcif_value(
                fields,
                indices,
                ["_atom_site.auth_seq_id", "_atom_site.label_seq_id"],
            ).strip()
            icode = _mmcif_value(fields, indices, ["_atom_site.pdbx_PDB_ins_code"]).strip()
            current_key = (chain_id, res_seq, icode)
            if current_key == last_key:
                continue
            last_key = current_key
            chains.setdefault(chain_id, []).append(three_to_one(resname))

    return chains


def parse_structure_chains(structure_path: Path) -> Dict[str, List[str]]:
    lower = structure_path.name.lower()
    if lower.endswith(".pdb"):
        return parse_pdb_chains(structure_path)
    if lower.endswith(".cif"):
        return parse_mmcif_chains(structure_path)
    raise ValueError(f"Unsupported structure format: {structure_path}")


def format_fasta(name: str, chains: Dict[str, List[str]]) -> str:
    entries = []
    multiple = len(chains) > 1

    for chain_id, residues in sorted(chains.items()):
        suffix = f"_{chain_id}" if multiple else ""
        header = f">{name}{suffix}"
        sequence = "".join(residues)
        entries.append(f"{header}\n{sequence}")

    return "\n".join(entries)


def write_fasta_file(output_dir: Path, pdb_id: str, chain_id: str, residues: List[str]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_chain = chain_id or "_"
    fasta_path = output_dir / f"{pdb_id}_{safe_chain}.fasta"
    header = f">{pdb_id}_{safe_chain}"
    sequence = "".join(residues)
    with fasta_path.open("w", encoding="ascii") as handle:
        handle.write(f"{header}\n{sequence}\n")
    return fasta_path


def extract_selected_chains(
    structure_path: Path,
    chain_ids: Sequence[str],
    output_dir: Path,
    label: str,
) -> None:
    try:
        chains = parse_structure_chains(structure_path)
    except ValueError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return

    if not chains:
        print(f"[warn] No residues parsed from {structure_path}", file=sys.stderr)
        return

    unique_chain_ids = []
    seen = set()
    for chain in chain_ids:
        normalized = (chain or "_").strip() or "_"
        if normalized not in seen:
            seen.add(normalized)
            unique_chain_ids.append(normalized)

    for chain_id in unique_chain_ids:
        residues = chains.get(chain_id)
        if residues is None:
            print(
                f"[warn] Chain {chain_id} not present in {structure_path}",
                file=sys.stderr,
            )
            continue
        output_path = write_fasta_file(output_dir, label, chain_id, residues)
        print(f"[info] Wrote {output_path}")


def run_default_batch() -> None:
    base = Path("pro")
    pdb_names = ["pro_a", "pro_b", "pro_c"]

    for pdb_name in pdb_names:
        pdb_file = base / f"{pdb_name}.pdb"
        if not pdb_file.exists():
            print(f"[warn] Missing file: {pdb_file}", file=sys.stderr)
            continue

        chains = parse_pdb_chains(pdb_file)
        if not chains:
            print(f"[warn] No chains found in {pdb_file}", file=sys.stderr)
            continue

        fasta_block = format_fasta(pdb_name, chains)
        print(fasta_block)


def main() -> None:
    args, parser = parse_args()

    if args.structure is not None:
        if not args.structure.exists():
            parser.error(f"Structure file not found: {args.structure}")
        if not args.chains:
            parser.error("At least one --chain value is required when --structure is set")
        label = args.label or args.structure.stem
        extract_selected_chains(args.structure, args.chains, args.output_dir, label)
        return

    if args.chains:
        parser.error("--chain requires --structure")

    run_default_batch()


if __name__ == "__main__":
    main()
