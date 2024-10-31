import sys
from lammps_methods import relax_box

if __name__ == '__main__':

    relax_box("fcc",3.615,10,"Cu_u3.eam",5000,63,180)
    sys.exit()

