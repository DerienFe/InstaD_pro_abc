#!/usr/bin/env python3
"""
Chain extractor for ProteinMPNN workflows.
- Selects specified chains from a PDB and writes a trimmed PDB.
- Emits chain sequences in order as a headerless '/'-joined string and (optionally) a small FASTA.
"""
import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils import seq1


def _clean_chain_ids(raw: str) -> List[str]:
    parts = raw.replace(',', ' ').split()
    uniq = []
    for p in parts:
        if p not in uniq:
            uniq.append(p)
    return uniq


def _chain_sequence(chain: PDB.Chain.Chain) -> str:
    # Build sequence directly from standard residues; include residues even if CA is missing so lengths match coords.
    seq_chars = []
    for res in chain:
        if res.id[0] != ' ':
            continue
        try:
            aa = seq1(res.resname)
        except Exception:
            continue
        if not aa or aa == 'X':
            continue
        seq_chars.append(aa)
    return ''.join(seq_chars)


def _count_standard_res(chain: PDB.Chain.Chain) -> int:
    return sum(1 for res in chain if res.id[0] == ' ')


def _residue_map(chain_ids: List[str], structure: PDB.Structure.Structure, source_path: Path) -> Dict[str, object]:
    """Build per-residue map (CA-only) with pdb_idx, chain_idx, global_idx, icode, aa."""
    rows: Dict[str, List[Dict[str, object]]] = {}
    global_idx = 1
    for cid in chain_ids:
        chain = structure[0][cid]
        chain_rows: List[Dict[str, object]] = []
        chain_idx = 1
        for res in chain:
            if res.id[0] != ' ':
                continue
            if 'CA' not in res:
                continue
            try:
                aa = seq1(res.resname)
            except Exception:
                continue
            if not aa or aa == 'X':
                continue
            pdb_idx = int(res.id[1])
            icode = res.id[2].strip() or ''
            chain_rows.append({
                'chain': cid,
                'global_idx': global_idx,
                'chain_idx': chain_idx,
                'pdb_idx': pdb_idx,
                'icode': icode,
                'aa': aa,
            })
            global_idx += 1
            chain_idx += 1
        rows[cid] = chain_rows
    return {
        'pdb': str(source_path),
        'scheme': 'pdb',
        'chains': chain_ids,
        'map': {cid: {'length': len(rlist), 'rows': rlist} for cid, rlist in rows.items()},
    }


def _write_map_files(map_payload: Dict[str, object], out_base: Path) -> None:
    base_no_ext = out_base.with_suffix('')
    out_json = base_no_ext.with_name(base_no_ext.name + '_map.json')
    out_csv = base_no_ext.with_name(base_no_ext.name + '_map.csv')
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(map_payload, indent=2))

    lines = ['chain,global_idx,chain_idx,pdb_idx,icode,aa']
    for cid in map_payload['chains']:
        for row in map_payload['map'][cid]['rows']:
            lines.append(f"{row['chain']},{row['global_idx']},{row['chain_idx']},{row['pdb_idx']},{row['icode']},{row['aa']}")
    out_csv.write_text('\n'.join(lines) + '\n')


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

        def accept_atom(self, atom):
            # Drop hetero atoms (including ligands, waters) to prevent sequence length mismatches.
            res = atom.get_parent()
            return res.id[0] == ' '

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

    # Emit residue map (CA-only) alongside outputs for downstream annotations.
    map_payload = _residue_map(chain_ids, structure, out_pdb)
    _write_map_files(map_payload, out_pdb)


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
