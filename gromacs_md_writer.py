#!/usr/bin/env python

def local_run_writer(workdir, mol_name):
    first_line = ("#!/bin/bash\n"
                  "export GMX_FORCE_UPDATE_DEFAULT_GPU=true\n")
    gromacs_run = ("gmx editconf -f complex.gro"
                   " -o box.gro -bt dodecahedron -c\n"
                   "gmx solvate -cp box.gro "
                   "-cs spc216.gro -o solv.gro -p topol.top\n"
                   "gmx grompp -f ions.mdp "
                   "-c solv.gro -p topol.top -o ions.tpr\n"
                   "echo SOL | "
                   "gmx genion -s ions.tpr -o solv_ions.gro"
                   " -p topol.top -pname NA "
                   "-nname CL -neutral -conc 0.15\n"
                   "# EM\n"
                   "gmx grompp -f em_steep.mdp "
                   "-c solv_ions.gro -p topol.top -o em.tpr -maxwarn 3\n"
                   "gmx mdrun -nb gpu -deffnm em -ntmpi 1 -ntomp 18 -pin on\n"
                   "# NVT\n"
                   "gmx grompp -f nvt.mdp "
                   "-c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 3\n"
                   "gmx mdrun -nb gpu -deffnm nvt -ntmpi 1 -ntomp 18 -pin on\n"
                   "# NPT\n"
                   "gmx grompp -f npt.mdp "
                   "-c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 3\n"
                   "gmx mdrun -deffnm npt -pin on -nb gpu -bonded gpu "
                   "-pme gpu -nstlist 400 -ntmpi 1 -ntomp 18\n"
                   "# Production MD\n"
                   "gmx grompp -f npt.mdp "
                   "-c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn 3\n"
                   "gmx mdrun -deffnm md -pin on -nb gpu -bonded gpu "
                   "-pme gpu -nstlist 400 -ntmpi 1 -ntomp 18\n")
    finish = ('echo -e "1 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 \n q" '
              '| gmx make_ndx -f md.gro -o index.ndx\n'
              'echo 1 0 | gmx trjconv -f md.xtc -o md_c.xtc -s md.tpr '
              '-pbc mol -center\n')
    jobruns = ("mpirun gmx_MMPBSA MPI -O -i mmgbsa.in "
               "-cs md.tpr -ci index.ndx"
               " -cg 42 22 -ct md_c.xtc "
               f"-cp topol.top -o {mol_name}.dat -nogui\n")

    with open(f"{workdir}/run_local.sh", "w") as f:
        f.write(first_line)
        f.write(gromacs_run)
        f.write(finish)
        f.write(jobruns)

    return 0

if __name__ == "__main__":
    local_run_writer('./', )