#this is a combined util file for all the functions used in the project
#by TW 10th Dec 2023

import os
import shutil
import subprocess
import config
import time
import sys
import argparse

import matplotlib.pyplot as plt

from glob import glob
from matplotlib import rcParams

rcParams.update({'font.size': 16})

def generate_tleap_scripts(parent_dir, gen_amber_tleap = False):

    # Define a template for the tleap script
    #dont change the indentation or enter new line as it perturbs the tleap script
    tleap_template = """source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p

pro_a = loadpdb ./pro/pro_a.amber.pdb
pro_b = loadpdb ./pro/pro_b.amber.pdb
pro_c = loadpdb ./pro/pro_c.amber.pdb

complex_ab = combine {pro_a pro_b}
complex_ac = combine {pro_a pro_c}

charge complex_ab
charge complex_ac

savepdb complex_ab ./mmgbsa/complex_ab.pdb
savepdb complex_ac ./mmgbsa/complex_ac.pdb


solvatebox complex_ab TIP3PBOX 12.0
solvatebox complex_ac TIP3PBOX 12.0

addions complex_ab Na+ 0
addions complex_ab Cl- 0
addions complex_ac Na+ 0
addions complex_ac Cl- 0

savepdb pro_a ./mmgbsa/ab/receptor.pdb
savepdb pro_a ./mmgbsa/ac/receptor.pdb
savepdb pro_b ./mmgbsa/ab/binder.pdb
savepdb pro_c ./mmgbsa/ac/binder.pdb

saveamberparm complex_ab ./mmgbsa/ab/complex.prmtop ./mmgbsa/ab/complex.inpcrd
saveamberparm complex_ac ./mmgbsa/ac/complex.prmtop ./mmgbsa/ac/complex.inpcrd
saveamberparm pro_a ./mmgbsa/ab/receptor.prmtop ./mmgbsa/ab/receptor.inpcrd
saveamberparm pro_a ./mmgbsa/ac/receptor.prmtop ./mmgbsa/ac/receptor.inpcrd
saveamberparm pro_b ./mmgbsa/ab/binder.prmtop ./mmgbsa/ab/binder.inpcrd
saveamberparm pro_c ./mmgbsa/ac/binder.prmtop ./mmgbsa/ac/binder.inpcrd
quit
"""
    #we write the template to output directory
    tleap_file_path = os.path.join(parent_dir, "gen_complex.tleap")
    with open(tleap_file_path, "w") as f:
        f.write(tleap_template)

    if gen_amber_tleap:
        raise NotImplementedError("amber tleap script generation is scrapped")

    def generate_minimization_input():
        minimization_template = """minimize complex
    &cntrl
    imin=1,
    maxcyc=1000,
    ncyc=500,
    cut=8.0,
    ntb=1,
    ntc=2,
    ntf=2,
    ntpr=100,
    ntr=1, 
    /"""
        return minimization_template

    def generate_NVT_equilibration_input():
        NVT_equilibration_template = """heat complex
    &cntrl
    imin=0,irest=0,ntx=1,
    nstlim=25000,dt=0.002,
    ntc=2,ntf=2,
    cut=8.0, ntb=1,
    ntpr=500, ntwx=500,
    ntt=3, gamma_ln=2.0,
    tempi=0.0, temp0=300.0, ig=-1,
    ntr=1, restraintmask=':1-242',
    restraint_wt=2.0,
    nmropt=1
    /
    &wt TYPE='TEMP0', istep1=0, istep2=25000,
    value1=0.1, value2=300.0, /
    &wt TYPE='END' /
    """
        return NVT_equilibration_template

    def generate_NPT_equilibration_input():
        NPT_equilibration_template = """heat complex
    &cntrl
    imin=0,irest=1,ntx=5,
    nstlim=25000,dt=0.002,
    ntc=2,ntf=2,
    cut=8.0, ntb=2, ntp=1, taup=1.0,
    ntpr=500, ntwx=500,
    ntt=3, gamma_ln=2.0,
    temp0=300.0, ig=-1,
    ntr=1, restraintmask=':1-242',
    restraint_wt=2.0,
    /"""
        return NPT_equilibration_template

    def generate_production_input():
        production_template = """heat complex
    &cntrl
    imin=0,irest=1,ntx=5,
    nstlim=250000,dt=0.002,
    ntc=2,ntf=2,
    cut=8.0, ntb=2, ntp=1, taup=2.0,
    ntpr=1000, ntwx=1000,
    ntt=3, gamma_ln=2.0,
    temp0=300.0, ig=-1,
    /"""
        return production_template

def convert_gmx(workdir, molname):
    import parmed as pmd
    # Convert molecule to Gromacs readible file
    prmtop_path = os.path.join(workdir, f'{molname}.prmtop')
    inpcrd_path = os.path.join(workdir, f'{molname}.inpcrd')
    par = pmd.load_file(prmtop_path, inpcrd_path)
    par.save(f'{workdir}/{molname}.gro', overwrite=True)
    par.save(f'{workdir}/{molname}.top', overwrite=True)
    return 0

def gmx_mdp_writer(workdir):
    ion_mdp = ('title		    = Ions\n'
               'integrator	    = steep	\n'
               'emtol		    = 1000.0  \n'
               'emstep          = 0.01     \n'
               'nsteps		    = 50000	  \n'
               'nstlist		    = 1		   \n'
               'cutoff-scheme   = Verlet\n'
               'ns_type		    = grid	\n'
               'rlist		    = 1.0	\n'
               'coulombtype	    = cutoff\n'
               'rcoulomb	    = 1.0	\n'
               'rvdw		    = 1.0	\n'
               'pbc             = xyz 	\n')
    em_mdp = ('define                 = -DFLEXIBLE\n'
              'integrator             = steep\n'
              'nsteps                 = 10000\n'
              'emtol                  = 100\n'
              'emstep                 = 0.01\n'
              'nstcomm                = 100\n'
              'nstxout                = 250        \n'
              'nstvout                = 0          \n'
              'nstfout                = 0          \n'
              'nstxout-compressed     = 100        \n'
              'compressed-x-precision = 1000\n'
              'nstlog                 = 500        \n'
              'nstenergy              = 500        \n'
              'nstcalcenergy          = 100\n'
              'cutoff-scheme          = Verlet\n'
              'ns-type                = grid\n'
              'nstlist                = 1\n'
              'rlist                  = 1.0\n'
              'constraints            = h-bonds\n'
              'coulombtype            = PME\n'
              'rcoulomb               = 1.0\n'
              'pme-order              = 6 \n'
              'fourierspacing         = 0.10\n'
              'ewald-rtol             = 1e-6\n'
              'vdw-type                = PME\n'
              'rvdw                    = 1.0\n'
              'vdw-modifier            = Potential-Shift\n'
              'ewald-rtol-lj           = 1e-3\n'
              'lj-pme-comb-rule        = Geometric\n'
              'DispCorr                = EnerPres\n'
              'Tcoupl              = no\n'
              'Pcoupl              = no\n'
              'gen_vel             = no\n')
    nvt_mdp = ('define       = -DPOSRE\n'
               'integrator   = sd      \n'
               'nsteps       = 50      \n'
               'dt           = 0.002   \n'
               'comm-mode    = Linear  \n'
               'nstcomm      = 100     \n'
               'nstxout                = 0      \n'
               'nstvout                = 0      \n'
               'nstfout                = 0      \n'
               'nstxout-compressed     = 5000   \n'
               'compressed-x-precision = 1000   \n'
               'nstlog                 = 5000   \n'
               'nstenergy              = 5000   \n'
               'nstcalcenergy          = 100    \n'
               'constraint_algorithm   = lincs  \n'
               'constraints            = h-bonds\n'
               'lincs_iter             = 1      \n'
               'lincs_order            = 4      \n'
               'lincs-warnangle        = 30     \n'
               'continuation           = no \n'
               'cutoff-scheme   = Verlet\n'
               'ns-type         = grid   \n'
               'nstlist         = 10     \n'
               'rlist           = 1.0    \n'
               'pbc             = xyz    \n'
               'coulombtype      = PME   \n'
               'rcoulomb         = 1.0   \n'
               'ewald_geometry   = 3d    \n'
               'pme-order        = 6     \n'
               'fourierspacing   = 0.10  \n'
               'ewald-rtol       = 1e-6  \n'
               'vdw-type                = PME\n'
               'rvdw                    = 1.0\n'
               'vdw-modifier            = Potential-Shift\n'
               'ewald-rtol-lj           = 1e-3\n'
               'lj-pme-comb-rule        = Geometric\n'
               'DispCorr                = EnerPres\n'
               'tc_grps          = System\n'
               'tau_t            = 1.0\n'
               'ref_t            = 300\n'
               'pcoupl           =  no\n'
               'gen_vel      = yes\n'
               'gen_seed     = -1 \n'
               'gen_temp     = 300\n')
    npt_mdp = ('define       = -DPOSRES\n'
               'integrator   = sd            \n'
               'nsteps       = 50000         \n'
               'dt           = 0.002         \n'
               'comm-mode    = Linear        \n'
               'nstcomm      = 100           \n'
               'nstxout                = 0       \n'
               'nstvout                = 0       \n'
               'nstfout                = 0       \n'
               'nstxout-compressed     = 100     \n'
               'compressed-x-precision = 1000    \n'
               'nstlog                 = 5000    \n'
               'nstenergy              = 500     \n'
               'nstcalcenergy          = 100     \n'
               'constraint_algorithm   = lincs   \n'
               'constraints            = h-bonds \n'
               'lincs_iter             = 1       \n'
               'lincs_order            = 4       \n'
               'lincs-warnangle        = 30      \n'
               'continuation           = yes     \n'
               'cutoff-scheme   = Verlet\n'
               'ns-type         = grid  \n'
               'nstlist         = 10    \n'
               'rlist           = 1.0   \n'
               'pbc             = xyz   \n'
               'coulombtype      = PME  \n'
               'rcoulomb         = 1.0  \n'
               'ewald_geometry   = 3d   \n'
               'pme-order        = 4    \n'
               'fourierspacing   = 0.10 \n'
               'ewald-rtol       = 1e-6 \n'
               'vdw-type                = cut-off\n'
               'rvdw                    = 1.0\n'
               'vdw-modifier            = Potential-Shift\n'
               'ewald-rtol-lj           = 1e-3\n'
               'DispCorr                = EnerPres\n'
               'tc_grps          = System\n'
               'tau_t            = 1.0\n'
               'ref_t            = 300\n'
               'pcoupl           = Berendsen\n'
               'pcoupltype       = isotropic\n'
               'tau_p            = 0.5                  \n'
               'ref_p            = 1.0                  \n'
               'compressibility  = 4.5e-05              \n'
               'refcoord-scaling = all\n'
               'gen_vel      = no\n')
    md_mdp = ('integrator   = sd   \n'
              'nsteps       = 500000        \n'
              'dt           = 0.002         \n'
              'comm-mode    = Linear        \n'
              'nstcomm      = 100           \n'
              'nstxout                = 0        \n'
              'nstvout                = 0        \n'
              'nstfout                = 0        \n'
              'nstxout-compressed     = 1000     \n'
              'compressed-x-precision = 1000     \n'
              'nstlog                 = 1000     \n'
              'nstenergy              = 1000     \n'
              'nstcalcenergy          = 100      \n'
              'constraint_algorithm   = lincs    \n'
              'constraints            = h-bonds  \n'
              'lincs_iter             = 1        \n'
              'lincs_order            = 4        \n'
              'lincs-warnangle        = 30       \n'
              'continuation           = yes      \n'
              'cutoff-scheme   = Verlet\n'
              'ns-type         = grid   \n'
              'nstlist         = 10     \n'
              'rlist           = 1.0    \n'
              'pbc             = xyz    \n'
              'coulombtype      = PME   \n'
              'rcoulomb         = 1.0   \n'
              'ewald_geometry   = 3d    \n'
              'pme-order        = 4     \n'
              'fourierspacing   = 0.10  \n'
              'ewald-rtol       = 1e-6  \n'
              'vdw-type                = cut-off\n'
              'rvdw                    = 1.0\n'
              'vdw-modifier            = Potential-Shift\n'
              'ewald-rtol-lj           = 1e-3\n'
              'DispCorr                = EnerPres\n'
              'tc_grps          = System\n'
              'tau_t            = 1.0\n'
              'ref_t            = 300\n'
              'pcoupl           = Parrinello-Rahma\n'
              'pcoupltype       = isotropic        \n'
              'tau_p            = 2                \n'
              'ref_p            = 1.0              \n'
              'compressibility  = 4.5e-05          \n'
              'gen_vel      = no\n')

    with open(f"{workdir}/ions.mdp", "w") as f:
        f.write(ion_mdp)
    with open(f"{workdir}/em_steep.mdp", "w") as f:
        f.write(em_mdp)
    with open(f"{workdir}/nvt.mdp", "w") as f:
        f.write(nvt_mdp)
    with open(f"{workdir}/npt.mdp", "w") as f:
        f.write(npt_mdp)
    with open(f"{workdir}/md.mdp", "w") as f:
        f.write(md_mdp)

    return 0

def mmgbsa_writer(workdir, type="gb"):
    qm = ('&general\n'
          'sys_name="QM-MMPBSA Decomposition",\n'
          'startframe=100,\n'
          'endframe=500,\n'
          'interval=20,\n'
          'forcefields="leaprc.protein.ff14SB,\n"'
          'temperature=300,\n'
          'interaction_entropy=1, ie_segment=25,\n'
          'c2_entropy=1\n'
          '/\n'
          '&gb\n'
          'igb=7, saltcon=0.150,\n'
          'ifqnt=1,\n'
          'qm_theory=PM6-DH+,\n'
          'qm_residues="within 3"\n'
          '/\n')
    pb_full = ('# General namelist variables \n'
            '&general \n'
            '  sys_name             = ""                                             # System name \n'
            '  startframe           = 100                                            # First frame to analyze \n'
            '  endframe             = 500                                            # Last frame to analyze \n'
            '  interval             = 20                                             # Number of frames between adjacent frames analyzed \n'
            '  forcefields          = "leaprc.ff14SB,leaprc.gaff"                    # Define the force field to build the Amber topology \n'
            '  ions_parameters      = 1                                              # Define ions parameters to build the Amber topology \n'
            '  PBRadii              = 3                                              # Define PBRadii to build amber topology from GROMACS files \n'
            '  temperature          = 298.15                                         # Temperature \n'
            '  qh_entropy           = 0                                              # Do quasi-harmonic calculation \n'
            '  interaction_entropy  = 0                                              # Do Interaction Entropy calculation \n'
            '  ie_segment           = 25                                             # Trajectory segment to calculate interaction entropy \n'
            '  c2_entropy           = 0                                              # Do C2 Entropy calculation \n'
            '  assign_chainID       = 0                                              # Assign chains ID \n'
            '  exp_ki               = 0.0                                            # Experimental Ki in nM \n'
            '  full_traj            = 0                                              # Print a full traj. AND the thread trajectories \n'
            '  gmx_path             = ""                                             # Force to use this path to get GROMACS executable \n'
            '  keep_files           = 2                                              # How many files to keep after successful completion \n'
            '  netcdf               = 0                                              # Use NetCDF intermediate trajectories \n'
            '  solvated_trajectory  = 1                                              # Define if it is necessary to cleanup the trajectories \n'
            '  verbose              = 1                                              # How many energy terms to print in the final output \n'
            '/ \n'
            '# (AMBER) Generalized-Born namelist variables \n'
            '&gb \n'
            '  igb                  = 5                                              # GB model to use \n'
            '  intdiel              = 1.0                                            # Internal dielectric constant for sander \n'
            '  extdiel              = 78.5                                           # External dielectric constant for sander \n'
            '  saltcon              = 0.150                                            # Salt concentration (M) \n'
            '  surften              = 0.0072                                         # Surface tension \n'
            '  surfoff              = 0.0                                            # Surface tension offset \n'
            '  molsurf              = 0                                              # Use Connelly surface (molsurf program) \n'
            '  msoffset             = 0.0                                            # Offset for molsurf calculation \n'
            '  probe                = 1.4                                            # Solvent probe radius for surface area calc \n'
            '  ifqnt                = 0                                              # Use QM on part of the system \n'
            '  qm_theory            = ""                                             # Semi-empirical QM theory to use \n'
            '  qm_residues          = ""                                             # Residues to treat with QM \n'
            '  qmcharge_com         = 0                                              # Charge of QM region in complex \n'
            '  qmcharge_lig         = 0                                              # Charge of QM region in ligand \n'
            '  qmcharge_rec         = 0                                              # Charge of QM region in receptor \n'
            '  qmcut                = 9999.0                                         # Cutoff in the QM region \n'
            '  scfconv              = 1e-08                                          # Convergence criteria for the SCF calculation, in kcal/mol \n'
            '  peptide_corr         = 0                                              # Apply MM correction to peptide linkages \n'
            '  writepdb             = 1                                              # Write a PDB file of the selected QM region \n'
            '  verbosity            = 0                                              # Controls the verbosity of QM/MM related output \n'
            '  alpb                 = 0                                              # Use Analytical Linearized Poisson-Boltzmann (ALPB) \n'
            '  arad_method          = 1                                              # Selected method to estimate the effective electrostatic size \n'
            '/ \n'
            '# (AMBER) Possion-Boltzmann namelist variables \n'
            '&pb \n'
            '  ipb                  = 2                                              # Dielectric model for PB \n'
            '  inp                  = 1                                              # Nonpolar solvation method \n'
            '  sander_apbs          = 0                                              # Use sander.APBS? \n'
            '  indi                 = 1.0                                            # Internal dielectric constant \n'
            '  exdi                 = 80.0                                           # External dielectric constant \n'
            '  emem                 = 4.0                                            # Membrane dielectric constant \n'
            '  smoothopt            = 1                                              # Set up dielectric values for finite-difference grid edges that are located across the solute/solvent dielectric boundary \n'
            '  istrng               = 0.0                                            # Ionic strength (M) \n'
            '  radiopt              = 1                                              # Use optimized radii? \n'
            '  prbrad               = 1.4                                            # Probe radius \n'
            '  iprob                = 2.0                                            # Mobile ion probe radius (Angstroms) for ion accessible surface used to define the Stern layer \n'
            '  sasopt               = 0                                              # Molecular surface in PB implict model \n'
            '  arcres               = 0.25                                           # The resolution (Å) to compute solvent accessible arcs \n'
            '  memopt               = 0                                              # Use PB optimization for membrane \n'
            '  mprob                = 2.7                                            # Membrane probe radius in Å \n'
            '  mthick               = 40.0                                           # Membrane thickness \n'
            '  mctrdz               = 0.0                                            # Distance to offset membrane in Z direction \n'
            '  poretype             = 1                                              # Use exclusion region for channel proteins \n'
            '  npbopt               = 0                                              # Use NonLinear PB solver? \n'
            '  solvopt              = 1                                              # Select iterative solver \n'
            '  accept               = 0.001                                          # Sets the iteration convergence criterion (relative to the initial residue) \n'
            '  linit                = 1000                                           # Number of SCF iterations \n'
            '  fillratio            = 4.0                                            # Ratio between the longest dimension of the rectangular finite-difference grid and that of the solute \n'
            '  scale                = 2.0                                            # 1/scale = grid spacing for the finite difference solver (default = 1/2 Å) \n'
            '  nbuffer              = 0.0                                            # Sets how far away (in grid units) the boundary of the finite difference grid is away from the solute surface \n'
            '  nfocus               = 2                                              # Electrostatic focusing calculation \n'
            '  fscale               = 8                                              # Set the ratio between the coarse and fine grid spacings in an electrostatic focussing calculation \n'
            '  npbgrid              = 1                                              # Sets how often the finite-difference grid is regenerated \n'
            '  bcopt                = 5                                              # Boundary condition option \n'
            '  eneopt               = 2                                              # Compute electrostatic energy and forces \n'
            '  frcopt               = 0                                              # Output for computing electrostatic forces \n'
            '  scalec               = 0                                              # Option to compute reaction field energy and forces \n'
            '  cutfd                = 5.0                                            # Cutoff for finite-difference interactions \n'
            '  cutnb                = 0.0                                            # Cutoff for nonbonded interations \n'
            '  nsnba                = 1                                              # Sets how often atom-based pairlist is generated \n'
            '  decompopt            = 2                                              # Option to select different decomposition schemes when INP = 2 \n'
            '  use_rmin             = 1                                              # The option to set up van der Waals radii \n'
            '  sprob                = 0.557                                          # Solvent probe radius for SASA used to compute the dispersion term \n'
            '  vprob                = 1.3                                            # Solvent probe radius for molecular volume (the volume enclosed by SASA) \n'
            '  rhow_effect          = 1.129                                          # Effective water density used in the non-polar dispersion term calculation \n'
            '  use_sav              = 1                                              # Use molecular volume (the volume enclosed by SASA) for cavity term calculation \n'
            '  cavity_surften       = 0.0378                                         # Surface tension \n'
            '  cavity_offset        = -0.5692                                        # Offset for nonpolar solvation calc \n'
            '  maxsph               = 400                                            # Approximate number of dots to represent the maximum atomic solvent accessible surface \n'
            '  maxarcdot            = 1500                                           # Number of dots used to store arc dots per atom \n'
            '  npbverb              = 0                                              # Option to turn on verbose mode \n'
            '/ \n'
            '# Decomposition namelist variables \n'
            '&decomposition \n'
            '  idecomp              = 3                                              # Which type of decomposition analysis to do \n'
            '  dec_verbose          = 3                                              # Control energy terms are printed to the output \n'
            '  print_res            = "within 6"                                     # Which residues to print decomposition data for \n'
            '  csv_format           = 1                                              # Write decomposition data in CSV format \n'
            '/ \n')
    gb = ('Input file for running PB and GB \n'
        'Taken from AMBER 2022 mannual page 853 \n'
        '&general \n'
        'sys_name="SARS_CoV2_S1 decomposition" \n'
        '  startframe=100, endframe=500, interval=20, \n'
        '  forcefields="leaprc.protein.ff14SB" \n'
        '  verbose=2, keep_files=0, \n'
        '/ \n'
        '&gb \n'
        '  igb=5, saltcon=0.150 \n'
        '/ \n'
        '&decomp \n'
        '  idecomp=1, dec_verbose=0, \n'
        '  print_res="within 4" \n'
        ' \n/')
    if type == "qm":
        with open(f"{workdir}/mmgbsa.in", "w") as f:
            f.write(qm)
    elif type == 'pb':
        with open(f"{workdir}/mmgbsa.in", "w") as f:
            f.write(pb_full)
    elif type == 'gb':
        with open(f"{workdir}/mmgbsa.in", "w") as f:
            f.write(gb)
    return 0

def local_run_writer(workdir, binder_name):
    """
    workdir: the directory where the files are located
    mol_name: the name of the molecule
    num_atom_pro_abc: the number of atoms in the protein A, B, C in list [A, B, C]
    """
    first_line = ("#!/bin/bash\n"
                  "export GMX_FORCE_UPDATE_DEFAULT_GPU=true\n")
    gromacs_run = ("# EM\n"
                   "gmx grompp -f em_steep.mdp "
                   "-c complex.gro -p complex.top -o em.tpr -maxwarn 3\n"
                   "gmx mdrun -nb gpu -deffnm em -ntmpi 1 -ntomp 18 -pin on\n"
                   "# NVT\n"
                   "gmx grompp -f nvt.mdp "
                   "-c em.gro -r em.gro -p complex.top -o nvt.tpr -maxwarn 3\n"
                   "gmx mdrun -nb gpu -deffnm nvt -ntmpi 1 -ntomp 18 -pin on\n"
                   "# NPT\n"
                   "gmx grompp -f npt.mdp "
                   "-c nvt.gro -r nvt.gro -p complex.top -o npt.tpr -maxwarn 3\n"
                   "gmx mdrun -deffnm npt -pin on -nb gpu -bonded gpu "
                   "-pme gpu -nstlist 400 -ntmpi 1 -ntomp 18\n"
                   "# Production MD\n"
                   "gmx grompp -f npt.mdp "
                   "-c npt.gro -r npt.gro -p complex.top -o md.tpr -maxwarn 3\n"
                   "gmx mdrun -deffnm md -pin on -nb gpu -bonded gpu "
                   "-pme gpu -nstlist 400 -ntmpi 1 -ntomp 18\n")
    finish = ('#grep the atom number use grep. \n'
              'atom_num_pro_receptor=$(grep -c "ATOM" receptor.pdb) \n'
              'atom_num_pro_binder=$(grep -c "ATOM" binder.pdb) \n'
              'start_pro_binder=$((atom_num_pro_receptor + 1)) \n'
              'end_pro_binder=$((atom_num_pro_receptor + atom_num_pro_binder)) \n'
              'echo -e "keep 1\\na 1-${atom_num_pro_receptor}\\nname 1 receptor\\na ${start_pro_binder}-${end_pro_binder}\\nname 2 binder\\nq" | gmx make_ndx -f md.gro -o index.ndx  \n'
              'echo 1 0 | gmx trjconv -f md.xtc -o md_c.xtc -s md.tpr '
              '-pbc mol -center\n')
    jobruns = ("#here we switch back to gmxMMPBSA conda env. \n"
               "source /home/tj/miniconda3/etc/profile.d/conda.sh  \n"
               "conda activate gmxMMPBSA  \n"
               "mpirun -np 16 gmx_MMPBSA MPI -O -i mmgbsa.in "
               "-cs md.tpr -ci index.ndx"
               " -cg 1 2 -ct md_c.xtc "
               f"-cp complex.top -nogui\n")

    with open(f"{workdir}/run_local.sh", "w") as f:
        f.write(first_line)
        f.write(gromacs_run)
        f.write(finish)
        f.write(jobruns)

    return 0

def file_preparation(mmgbsa_dir):
    """
    mmgbsa_dir: the directory where the files are located
    receptor: the name of the receptor [pro_a]
    binder: the name of the binder [pro_b, pro_c]
    """

    workdirs = glob(f"{mmgbsa_dir}/*")
    workdirs = [workdir for workdir in workdirs if os.path.isdir(workdir)]
    for workdir in workdirs:
        #conversion of prmtop/inpcrd to gro/top
        prmtop_files = glob(f"{workdir}/complex.prmtop")
        basenames = [os.path.basename(prmtop_file).split(".")[0] for prmtop_file in prmtop_files]
        for basename in basenames:
            convert_gmx(workdir, basename)
            print(f"converting {basename} to gmx format")
        
        #generate mdp files for running MD gromacs
        gmx_mdp_writer(workdir)
        print(f"generating mdp files for {workdir}...")

        #generate mmgbsa input files for running mmgbsa gromacs
        mmgbsa_writer(workdir)
        print(f"generating mmgbsa input files for {workdir}...")
    
        #generate local .sh file to run gromacs on the generated mdp files
        local_run_writer(workdir, basename)
    return 0

def read_decomp_results(datfile, exp_name='decomp_exp', receptor_offset=333, binder_offset=1, plot = True):
    with open(datfile,"r") as file:
        read_flag = False
        result_lines = []
        for line in file:
            if "Total Energy Decomposition:" in line:
                read_flag = True
            if read_flag:
                result_lines.append(line)
    
    receptor_lines = [line for line in result_lines if "R:A" in line]
    binder_lines = [line for line in result_lines if "L:B" in line]

    #we get the 2,3 element combined as string for label, and -3 element for value
    receptor_labels = [line.split(",")[0] for line in receptor_lines]
    receptor_values = [float(line.split(",")[-3]) for line in receptor_lines]
    binder_labels = [line.split(",")[0] for line in binder_lines]
    binder_values = [float(line.split(",")[-3]) for line in binder_lines]

    def fix_index(labels, offset):
        new_labels = []
        for label in labels:
            index = int(label.split(":")[-1])
            pre = label.split(":")[:-1]
            index = index + offset -1
            new_label = pre[-1] + ":" + str(index) #":".join(pre) + ":" + str(index)
            new_labels.append(new_label)
        return new_labels
    
    receptor_labels = fix_index(receptor_labels, receptor_offset)
    binder_labels = fix_index(binder_labels, binder_offset)

    if plot:
        figs_dir = os.path.join("./","fig")
        if not os.path.exists(figs_dir):
            os.makedirs(figs_dir)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        fig.suptitle(exp_name)
        ax1.bar(receptor_labels, receptor_values)
        ax1.set_title("Receptor decomp result")
        ax2.bar(binder_labels, binder_values)
        ax2.set_title("Binder decomp result")
        ax1.set_ylabel("kcal/mol")
        ax2.set_ylabel("kcal/mol")
        plt.setp(ax1.get_xticklabels(), rotation=60, fontsize=12)
        plt.setp(ax2.get_xticklabels(), rotation=60, fontsize=12)
        fig.savefig(os.path.join(figs_dir, f"{exp_name}.png"))
        plt.close(fig)

    return 0