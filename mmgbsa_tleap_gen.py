#generation of tleap scripts for mmgbsa
#by TW 07th Dec 2023

from util import *

if __name__ == "__main__":
    generate_tleap_scripts('./mmgbsa/')
    print("tleap scripts generated at time: ", time.ctime())

