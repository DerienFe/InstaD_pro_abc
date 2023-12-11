#this is the main sh file for whole automated process
#by TW 07th Dec 2023
#initialize all the folder structure

########################
# 0. initialization
########################

#!/bin/bash
mkdir ./mmgbsa
mkdir ./mmgbsa/ab
mkdir ./mmgbsa/ac

mkdir ./metaD
mkdir ./metaD/trajectory
mkdir ./metaD/aux_file_dir

########################
# 1. system preparation
########################

source /home/tj/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

pdb4amber -i ./pro/pro_a.pdb -o ./pro/pro_a.amber.pdb -y >> log.txt &&
pdb4amber -i ./pro/pro_b.pdb -o ./pro/pro_b.amber.pdb -y >> log.txt &&
pdb4amber -i ./pro/pro_c.pdb -o ./pro/pro_c.amber.pdb -y >> log.txt &&

python mmgbsa_tleap_gen.py >> log.txt &&

tleap -f ./mmgbsa/gen_complex.tleap >> log.txt &&

#convert prmtop/inpcrd to gro since we only have access to gromacs...
python file_preparation.py >> log.txt &&

########################
# 2. run simulation
########################

#deactivate gmxMMPBSA conda env. source the gromacs
#this is due to gmxMMPBSA was not complied with CUDA.
conda deactivate
source /usr/local/gromacs/bin/GMXRC

cd ./mmgbsa/ac
chmod +x run_local.sh
bash run_local.sh >> ../../log.txt &&

cd ../ab
chmod +x run_local.sh 
bash run_local.sh >> ../../log.txt &&

cd ../../

########################
# 3. post-MD analysis
########################

python mmgbsa_decomp_result_analysis.py >> log.txt

#result saved in fig/

#meta-Dynamics
conda deactivate
conda activate biophys_env

python metaD_sim.py >> log.txt

#result saved in metaD/