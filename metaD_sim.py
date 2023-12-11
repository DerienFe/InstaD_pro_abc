#metaDynamics simulation of Phe486 of Spike protein S1 SARS-CoV19 with Tyr52 on Ab P5A3C8 (THSC20.HVTR26) 

import numpy as np
import matplotlib.pyplot as plt
import os
import time

from tqdm import tqdm

import openmm as omm
from openmm import unit
from openmm import Vec3
import openmm.app as omm_app
from openmm.app import *
from openmm.app.topology import Topology
from openmm.app.element import Element
from openmm.app.metadynamics import BiasVariable, Metadynamics
from openmm.unit import Quantity

from util import *
import config_openmm as config
import csv

def minimize(context):
    st = time.time()
    s = time.time()
    print("Setting up the simulation")

    # Minimizing step
    context.setPositions(gro.positions)
    state = context.getState(getEnergy = True)
    energy = state.getPotentialEnergy()

    for _ in range(50):
        omm.openmm.LocalEnergyMinimizer.minimize(context, 1, 20)
        state = context.getState(getEnergy = True)
        energy = state.getPotentialEnergy()

    print("Minimization done in", time.time() - s, "seconds")
    s = time.time()
    return context, energy



if __name__ == "__main__":
    gro_file = './mmgbsa/ab/complex.gro'
    top_file = './mmgbsa/ab/complex.top'

    meta_freq = 5000
    meta_height = 5     #kcal/mol

    gro = GromacsGroFile(gro_file)
    top = GromacsTopFile(top_file, periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir='/usr/local/gromacs/share/gromacs/top')
    #forcefield = omm_app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    #params = omm_app.CharmmParameterSet('toppar/toppar_water_ions.str') #we modified the combined LJ term between NaCl to have a -6.0kcal.mol at 2.5A

    #save the pdb file for check.
    #with open('./metaD/complex.pdb', 'w') as f:
    #    omm_app.PDBFile.writeFile(prmtop.topology, inpcrd.positions, f)

    for i_sim in range(config.num_sim):
        print(f"MetaD Simulation {i_sim} starting")
        time_tag = time.strftime("%Y%m%d-%H%M%S")
        system = top.createSystem(nonbondedMethod=PME,
                                  nonbondedCutoff=1.0*unit.nanometers,
                                  constraints=HBonds)
        
        #metaD
        centroid_bond_force = omm.CustomCentroidBondForce(2, "distance(g1,g2)")
        group1_index = []
        group2_index = []
        for atom in top.topology.atoms():
            if atom.residue.index == config.res1_index:
                group1_index.append(atom.index)
            if atom.residue.index == config.res2_index:
                group2_index.append(atom.index)
        print("group1_index", group1_index)
        print("group2_index", group2_index)

        group1_index = centroid_bond_force.addGroup(group1_index)
        group2_index = centroid_bond_force.addGroup(group2_index)
        centroid_bond_force.addBond([group1_index, group2_index])
        dist_cv = BiasVariable(centroid_bond_force, 0.3, 1, 0.02, True)  #at PDB 7z0x distance ~= 4.8A. note openmm use nm.

        aux_file_path = os.path.join("./metaD/aux_file_dir/" + time_tag)
        if not os.path.exists(aux_file_path):
            os.makedirs(aux_file_path, exist_ok=True)
        metaD = Metadynamics(system=system, 
                            variables=[dist_cv], #define the cv
                            temperature=config.T*unit.kelvin,
                            biasFactor=5.0,
                            height=meta_height * 4.184 * unit.kilojoules_per_mole,
                            frequency=meta_freq,
                            saveFrequency=meta_freq,
                            biasDir=aux_file_path)
        
        platform = config.platform
        integrator = omm.LangevinIntegrator(config.T*unit.kelvin, #Desired Integrator
                                            10/unit.picoseconds,
                                            config.stepsize)

        #simulation object, note the context object is automatically created.
        sim = omm.app.Simulation(top.topology, system, integrator, platform=platform)
        sim.context.setPositions(gro.positions)

        #minimize the energy
        context, energy = minimize(sim.context)
        sim.context = context

        #run the simulation
        file_handle = open(f'./metaD/trajectory/{time_tag}_metaD_traj.dcd', 'wb')
        dcd_file = omm_app.DCDFile(file_handle, top.topology, dt = config.stepsize)

        for _ in tqdm(range(int(config.sim_steps/config.dcdfreq))):
            metaD.step(sim, config.dcdfreq)
            state = context.getState(getPositions=True)
            dcd_file.writeModel(state.getPositions(asNumpy=True))
        file_handle.close()
        energy_CV = metaD.getFreeEnergy() #The values are in kJ/mole. The i’th position along an axis corresponds to minValue + i*(maxValue-minValue)/gridWidth.
        gridwidth = int(np.ceil(5*(1.0-0.2)/0.02)) #200.
        cv_space = np.linspace(2, 10, gridwidth)

        #create the fes metad folder is not exist.
        if not os.path.exists("fes_metaD"):
            os.makedirs("fes_metaD")

        #plot the free energy surface
        plt.figure()
        plt.plot(cv_space, energy_CV/4.184)
        plt.xlabel("distance between COM of Phe486 and Tyr52")
        plt.ylabel("Free energy (kcal/mol)")
        plt.savefig(f"./metaD/{time_tag}_fes_metaD_{i_sim}.png")
        plt.close()

    
    
    
    
    
    
    
