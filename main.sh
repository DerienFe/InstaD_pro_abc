#this is the main sh file for whole automated process
#by TW 07th Dec 2023
#initialize all the folder structure

#!/bin/bash
mkdir ./prep
mkdir ./mmgbsa

mkdir ./metaD
mkdir ./metaD/trajectory
mkdir ./metaD/aux_file_dir

#we generate the tleap input file for mmgbsa
#python mmgbsa_tleap_gen.py 
#adapting Denes code to run EM, NVT, NPT, prod on GROMACS
#and generate 

#set pro folder name.

gmx pdb2gmx -f ./pro/pro_a.pdb -o ./prep/pro_a_processed.gro -p ./prep/pro_a.top -water spce
gmx pdb2gmx -f ./pro/pro_b.pdb -o ./prep/pro_b_processed.gro -p ./prep/pro_b.top -water spce
gmx pdb2gmx -f ./pro/pro_c.pdb -o ./prep/pro_c_processed.gro -p ./prep/pro_c.top -water spce

#iteratively get in the folder and call the input_writer and gromacs_md_writer

