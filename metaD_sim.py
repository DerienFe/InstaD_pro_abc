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
    context.setPositions(pdb.positions)
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
    prmtop_file = './mmgbsa/complex.prmtop' # Path #tpr #prmtop
    inpcrd_file = './mmgbsa/complex.inpcrd' # Path # gro # inpcrd

    meta_freq = 5000
    meta_height = 5     #kcal/mol

    prmtop = omm_app.AmberPrmtopFile(prmtop_file)
    inpcrd = omm_app.AmberInpcrdFile(inpcrd_file)
    #forcefield = omm_app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    params = omm_app.CharmmParameterSet('toppar/toppar_water_ions.str') #we modified the combined LJ term between NaCl to have a -6.0kcal.mol at 2.5A

    for i_sim in range(config.num_sim):
        print(f"MetaD Simulation {i_sim} starting")
        time_tag = time.strftime("%Y%m%d-%H%M%S")
        system = prmtop.createSystem(params,
                                  nonbondedCutoff=1.0*unit.nanometers,
                                  constraints=omm_app.HBonds)
        
        #metaD
        bias_bond = omm.CustomBondForce("r")
        bias_bond.addBond(0, 1)
        dist_cv = BiasVariable(bias_bond, 0.3, 1, 0.02, False)  #at PDB 7z0x distance = 4.8A. note openmm use nm.

        bias_dist = omm.CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)^2)")
        bias_dist.addGroup([486]) #index of Phe486
        bias_dist.addGroup([52]) #index of Tyr52 
        


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
        
        platform = omm.Platform.getPlatformByName('CUDA')
        integrator = omm.LangevinIntegrator(config.T*unit.kelvin, #Desired Integrator
                                            10/unit.picoseconds,
                                            config.stepsize)

        NaCl_dist = [[]] #initialise the NaCl distance list.
        #simulation object, note the context object is automatically created.
        sim = omm.app.Simulation(inpcrd.topology, system, integrator, platform=platform)
        sim.context.setPositions(inpcrd.positions)

        #minimize the energy
        context, energy = minimize(sim.context)
        sim.context = context

        #run the simulation
        file_handle = open(f'/metaD/trajectory/{time_tag}_metaD_traj.dcd', 'wb')
        dcd_file = omm_app.DCDFile(file_handle, inpcrd.topology, dt = config.stepsize)

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

    
    
    
    
    
    
    
