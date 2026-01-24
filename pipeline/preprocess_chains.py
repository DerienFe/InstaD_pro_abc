#!/usr/bin/env python3
"""
Chain extractor for ProteinMPNN workflows.
- Selects specified chains from a PDB and writes a trimmed PDB.
- Emits chain sequences in order as a headerless '/'-joined string and (optionally) a small FASTA.
"""
import argparse
from pathlib import Path
from typing import List

from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder


def _clean_chain_ids(raw: str) -> List[str]:
    parts = raw.replace(',', ' ').split()
    uniq = []
    for p in parts:
        if p not in uniq:
            uniq.append(p)
    return uniq


def _chain_sequence(chain: PDB.Chain.Chain) -> str:
    builder = PPBuilder()
    peptides = builder.build_peptides(chain)
    seq = ''.join(str(pp.get_sequence()) for pp in peptides)
    return seq


def _count_standard_res(chain: PDB.Chain.Chain) -> int:
    return sum(1 for res in chain if res.id[0] == ' ')


def extract_chains(pdb_path: Path, chain_ids: List[str], out_pdb: Path, out_fasta: Path | None, out_seq: Path | None) -> None:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(out_pdb.stem, pdb_path)
    available = {chain.id for chain in structure.get_chains()}
    missing = [c for c in chain_ids if c not in available]
    if missing:
        raise ValueError(f"Chains not found in {pdb_path.name}: {missing}")

    class ChainSelect(PDB.Select):
        def accept_chain(self, chain):
            return chain.id in chain_ids

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), ChainSelect())

    sequences = []
    for cid in chain_ids:
        chain = structure[0][cid]
        seq = _chain_sequence(chain)
        observed = _count_standard_res(chain)
        if observed != len(seq):
            raise ValueError(f"Chain {cid} residue count mismatch: coords {observed}, seq {len(seq)}")
        sequences.append(seq)

    joined = '/'.join(sequences)
    if out_seq:
        out_seq.write_text(joined + "\n")

    if out_fasta:
        header = f">{out_pdb.stem}|{''.join(chain_ids)}\n"
        out_fasta.write_text(header + joined + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Trim PDB to selected chains and emit sequences")
    parser.add_argument('--pdb', required=True, type=Path)
    parser.add_argument('--chains', required=True, help="Chain IDs to keep, space or comma separated")
    parser.add_argument('--out-pdb', required=True, type=Path)
    parser.add_argument('--out-fasta', type=Path, default=None)
    parser.add_argument('--out-seq', type=Path, default=None, help="Headerless '/'-joined sequence output")
    args = parser.parse_args()

    chain_ids = _clean_chain_ids(args.chains)
    extract_chains(args.pdb, chain_ids, args.out_pdb, args.out_fasta, args.out_seq)


if __name__ == '__main__':
    main()
