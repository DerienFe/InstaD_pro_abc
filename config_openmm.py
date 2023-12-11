import time
from openmm import unit
import openmm

from openmm.unit import Quantity
from openmm import Vec3

num_sim = 2
sim_steps = int(5e7)
pbc = False
time_tag = time.strftime("%Y%m%d-%H%M%S")
T = 300
stepsize = 0.002 * unit.picoseconds #equivalent to 2 * unit.femtoseconds 4fs.
dcdfreq = 5000
platform = openmm.Platform.getPlatformByName('CUDA')

res1_index = int(486-333) # Phe486, minus the 333 offset.
res2_index = int(52+194) # the S1 protein has total of 194 residues, plus the 52 residue number