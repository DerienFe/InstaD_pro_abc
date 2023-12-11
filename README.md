# Automated Pipeline for MMGBSA Decomposition and 1D Meta-Dynamics

## Introduction
This repository contains an automated pipeline designed for Molecular Mechanics Generalized Born Surface Area (MMGBSA) per-residue decomposition and 1D meta-Dynamics simulations. The pipeline is optimized for performance and accuracy, providing a comprehensive solution for researchers in the field of computational biophysics.

## Prerequisites
Before setting up the pipeline, ensure that you have the following prerequisites:
- A Linux operating system or a compatible bash shell environment.
- Gromacs with CUDA compilation support for enhanced performance.
- Conda package manager with two specific environments: `gmxMMPBSA` and `biophys_env`.

## Installation Guide
### Gromacs Installation
For detailed instructions on installing Gromacs, please refer to the [Gromacs Installation Guide](https://manual.gromacs.org/documentation/current/install-guide/index.html).
To install Gromacs with GPU support, follow these steps:
tar xfz gromacs-2023.3.tar.gz
cd gromacs-2023.3
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DGMX_MPI=on
make
make check
sudo make install
source /usr/local/gromacs/bin/GMXRC


### Setting Up Conda Environments
Conda environments are used to manage the dependencies and packages required for this pipeline. For instructions on installing Conda, visit the [Conda Installation Guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

#### gmxMMPBSA Environment
To set up the `gmxMMPBSA` environment, follow the guidelines provided [here](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/installation/).

#### biophys_env Environment
To create and configure the `biophys_env` environment, execute the following commands:
conda create -n biophys_env
conda activate biophys_env
conda install -c conda-forge openmm python numpy scipy tqdm scipy matplotlib mdtraj plotly python-kaleido


## Usage
chmod +x ./main.sh
bash main.sh

