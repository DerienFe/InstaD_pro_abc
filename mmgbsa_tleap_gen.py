#generation of tleap scripts for mmgbsa
#by TW 07th Dec 2023

import argparse
import os
import config

def generate_tleap_scripts(combinations, gen_tleap = True):
    # Create argument parser
    parser = argparse.ArgumentParser(description='Generate tleap scripts for mmgbsa')
    parser.add_argument('-od', '--output_directory', default='./mmgbsa', help='Output directory for tleap scripts')

    # Parse command line arguments
    args = parser.parse_args()

    # Define a template for the tleap script
    #dont change the indentation or enter new line as it perturbs the tleap script
    tleap_template = """{complex} = loadpdb ../pro/{complex_filename}.pdb
{component1} = loadpdb ../pro/{component1_filename}.pdb
{component2} = loadpdb ../pro/{component2_filename}.pdb
set default PBRadii mbondi2
saveamberparm {complex} {complex_filename}.prmtop {complex_filename}.inpcrd
saveamberparm {component1} {component1_filename}.prmtop {component1_filename}.inpcrd
saveamberparm {component2} {component2_filename}.prmtop {component2_filename}.inpcrd

charge {complex}
solvatebox {complex} TIP3Ppro_bOX 12.0
saveamberparm {complex} {complex_filename}_solvated.prmtop {complex_filename}_solvated.inpcrd
quit
"""

    # Define the different combinations
    

    # Generate tleap scripts for each combination
    tleap_scripts = {}

    for combo in combinations:
        complex_name, comp1, comp2 = combo
        tleap_script = tleap_template.format(
            complex=complex_name.replace("-", ""),
            complex_filename=complex_name,
            component1=comp1,
            component2=comp2,
            component1_filename=comp1,
            component2_filename=comp2
        )
        tleap_scripts[complex_name] = tleap_script


    # Create the subfolder for each complex_name
    for complex_name in tleap_scripts.keys():
        subfolder_path = os.path.join(args.output_directory, complex_name)
        os.makedirs(subfolder_path, exist_ok=True)

    if gen_tleap:
    # Write the tleap scripts to files within the subfolders
        for complex_name, tleap_script in tleap_scripts.items():
            # Define the subfolder path for each complex
            subfolder_path = os.path.join(args.output_directory, complex_name)

            # Write tleap script
            tleap_file_path = os.path.join(subfolder_path, f"{complex_name}.tleap")
            with open(tleap_file_path, "w") as f:
                f.write(tleap_script)
            
            # Generate and write minimization input
            min_file_path = os.path.join(subfolder_path, f"{complex_name}_min.in")
            with open(min_file_path, "w") as f:
                f.write(generate_minimization_input())

            # Generate and write NVT equilibration input
            nvt_file_path = os.path.join(subfolder_path, f"{complex_name}_NVT.in")
            with open(nvt_file_path, "w") as f:
                f.write(generate_NVT_equilibration_input())
            
            # Generate and write NPT equilibration input
            npt_file_path = os.path.join(subfolder_path, f"{complex_name}_NPT.in")
            with open(npt_file_path, "w") as f:
                f.write(generate_NPT_equilibration_input())
            
            # Generate and write production input
            production_file_path = os.path.join(subfolder_path, f"{complex_name}_production.in")
            with open(production_file_path, "w") as f:
                f.write(generate_production_input())

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


if __name__ == "__main__":
    combinations = config.combinations
    generate_tleap_scripts(combinations)
