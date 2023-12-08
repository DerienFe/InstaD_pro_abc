#Auto align function for A-B, A-C alignment using pymol API
#by Tiejun Wei 07th Dec 2023

from pymol import cmd
import os
import sys

def auto_align(pdb1, pdb2, chain1, chain2, binder1, binder2):
    """
    Usage: auto_align(pdb1, pdb2, chain1, chain2)
    pdb1: str path to the pdb file1
    pdb2: str path to the pdb file2
    chain1: str aligning chain name of the pdb1 (pro_A)
    chain2: str aligning chain name of the pdb2 (pro_A)
    binder1: str chain name of the binded ligand/receptor to pro_A in pdb1 (pro_B)
    binder2: str chain name of the binded ligand/receptor to pro_A in pdb2 (pro_C)
    note the pdb1 is moving, pdb2 is fixed
    please use the auth chain name specified in the PDB website
    """
    
    for pdb in [pdb1, pdb2]:
        assert os.path.exists(pdb), f"{pdb} does not exist, please check the path"
    
    #load the pdb files
    cmd.load(pdb1, "pdb1")
    cmd.load(pdb2, "pdb2")

    #we process the list and make it a selection string for pymol
    pro_A_chain_mask_pdb1 = " or ".join([f"chain {i.upper()}" for i in chain1])
    pro_A_chain_mask_pdb2 = " or ".join([f"chain {i.upper()}" for i in chain2])
    pro_B_chain_mask = " or ".join([f"chain {i.upper()}" for i in binder1])
    pro_C_chain_mask = " or ".join([f"chain {i.upper()}" for i in binder2])

    #select the chains
    cmd.select("chain1", f"pdb1 and ({pro_A_chain_mask_pdb1})")
    cmd.select("chain2", f"pdb2 and ({pro_A_chain_mask_pdb2})")
    cmd.select("binder1", f"pdb1 and ({pro_B_chain_mask})")
    cmd.select("binder2", f"pdb2 and ({pro_C_chain_mask})")

    #check if the select is empty
    for chain in ["chain1", "chain2"]:
        assert cmd.count_atoms(chain) != 0, f"{chain} is empty, please check the chain name"
    
    #align the chain
    cmd.align("chain1", "chain2")
    
    #merge the two pdb and save it.
    cmd.create("merged", "chain1 or chain2")
    cmd.save(f"./pro/{pdb1.split('/')[-1].split('.')[0]}_{pdb2.split('/')[-1].split('.')[0]}_merged.pdb", "merged")

    #also save the aligned pdb seperately into 3 parts: pro_a, pro_b, pro_c
    cmd.save(f"./pro/pro_a.pdb", "chain1")
    cmd.save(f"./pro/pro_b.pdb", "binder1")
    cmd.save(f"./pro/pro_c.pdb", "binder2")

    cmd.delete("all")
    cmd.quit()

if __name__ == "__main__":
    pdb1 = "./pro/6m0j.pdb"
    pdb2 = "./pro/7z0x.pdb"
    pro_A_chain_in_pdb1 = ["E"]
    pro_A_chain_in_pdb2 = ["R"]

    pro_B_chain_in_pdb1 = ["A"]
    pro_C_chain_in_pdb2 = ["H", "L"]

    auto_align(pdb1, pdb2, pro_A_chain_in_pdb1, pro_A_chain_in_pdb2, pro_B_chain_in_pdb1, pro_C_chain_in_pdb2)
    print("Done")
