import time
from openmm import unit
import openmm

from openmm.unit import Quantity
from openmm import Vec3

num_sim = 20
sim_steps = int(5e7)
pbc = False
time_tag = time.strftime("%Y%m%d-%H%M%S")
T = 300
stepsize = 0.002 * unit.picoseconds #equivalent to 2 * unit.femtoseconds 4fs.
dcdfreq = 500
platform = openmm.Platform.getPlatformByName('CUDA')