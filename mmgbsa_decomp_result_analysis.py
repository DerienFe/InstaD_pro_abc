#load the data and post-process it.
from util import *

if __name__ == "__main__":
    read_decomp_results('./mmgbsa/ab/FINAL_DECOMP_MMPBSA.dat', "S1 - Ab", receptor_offset=333, binder_offset=1,)
    read_decomp_results('./mmgbsa/ac/FINAL_DECOMP_MMPBSA.dat', "S1 - ACE2", receptor_offset=333, binder_offset=19)