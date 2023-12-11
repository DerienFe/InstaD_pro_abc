
# Automated Molecular Simulation Process
_By TW, 7th Dec 2023_

This document outlines the `main.sh` script used for automating the molecular simulation process.

## 0. Initialization
Initialize all the folder structure necessary for the process.

```bash
#!/bin/bash
mkdir ./mmgbsa
mkdir ./mmgbsa/ab
mkdir ./mmgbsa/ac

mkdir ./metaD
mkdir ./metaD/trajectory
mkdir ./metaD/aux_file_dir
```

## 1. System Preparation
Prepare the system by activating environments and processing PDB files.

```bash
source /home/tj/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

pdb4amber -i ./pro/pro_a.pdb -o ./pro/pro_a.amber.pdb -y >> log.txt &&
pdb4amber -i ./pro/pro_b.pdb -o ./pro/pro_b.amber.pdb -y >> log.txt &&
pdb4amber -i ./pro/pro_c.pdb -o ./pro/pro_c.amber.pdb -y >> log.txt &&

python mmgbsa_tleap_gen.py >> log.txt &&

tleap -f ./mmgbsa/gen_complex.tleap >> log.txt &&

# convert prmtop/inpcrd to gro since we only have access to gromacs...
python file_preparation.py >> log.txt &&
```

## 2. Run Simulation
Execute the simulation process by changing environments and running local scripts.

```bash
# deactivate gmxMMPBSA conda env. source the gromacs
# this is due to gmxMMPBSA was not complied with CUDA.
conda deactivate
source /usr/local/gromacs/bin/GMXRC

cd ./mmgbsa/ac
chmod +x run_local.sh
bash run_local.sh >> ../../log.txt &&

cd ../ab
chmod +x run_local.sh 
bash run_local.sh >> ../../log.txt &&

cd ../../
```

## 3. Post-MD Analysis
Conduct post-Molecular Dynamics analysis and metadynamics simulation.

```bash
python mmgbsa_decomp_result_analysis.py >> log.txt

# result saved in fig/

# meta-Dynamics
conda deactivate
conda activate biophys_env

python metaD_sim.py >> log.txt

# result saved in metaD/
```
