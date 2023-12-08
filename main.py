#this is the main python files taking over the original main.sh bash script
#by Tiejun Wei 08th Dec 2023

from input_writer import *
from gromacs_md_writer import *
from config import *
from glob import glob
import os
import shutil

def setup_parent_dir(dir):
    path = os.path.join(dir, 'tmpdir')
    os.makedirs(path, exist_ok=True)
    return path

def main(parent_dir,
         binder_mol,
         receptor_mol,
         output_dir,
         mode,
         parallel,
         num_cores,
         use_gpu,
         ):
    """
    parent_dir: the parent directory for all the calculations
    binder_mol: the path to the binder molecule
    receptor_mol: the path to the receptor molecule
    output_dir: the output directory
    mode: the calculation mode (PB, GB, QM)
    parallel: whether to use parallel or not
    num_cores: number of cores to use
    use_gpu: flag to use gpu or not
    """
    return 0
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='MMPBSA decomposition')
    parser.add_argument('-b', '--binder', help='list of binder protein (pro_b, pro_c)', default=['pro_b', 'pro_c'])
    parser.add_argument('-r', '--receptor', help='the receptor protein (pro_a)', default='pro_a')
    parser.add_argument('-id', '--input_dir', help='input directory', default='./pro/')
    parser.add_argument('-od', '--output_dir', help='output directory', default='./mmgbsa/')
    parser.add_argument('-m', '--mode', help='select from "PB", "GB" or "QM"', default='PB')
    parser.add_argument('-p', '--parallel', help='use parallel or not', default=False)
    parser.add_argument('-n', '--num', help='number of cores', default=8)
    parser.add_argument('-g', '--gpu', help='use gpu or not', default=False)
    parse = parser.parse_args()

    parent_dir = setup_parent_dir(parse.output_dir)
    
    binder_mol =  os.path.join(parse.input_dir, parse.binder)
    receptor_mol = os.path.join(parse.input_dir, parse.receptor)

    main()

    print("preparation finished")