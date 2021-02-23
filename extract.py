#!/bin/env python3


from dispersion import *
import sys
import os

def main():

    path   = sys.argv[1]
    ncells = int(sys.argv[2])
    component = sys.argv[3]
    t = get_times(os.path.join(path, "EM_B.h5"))
    r = Run(path)
    B = mergeAllTimes(r, ncells, t, component=component, silent=False)
    save(B, os.path.join(path,"{}t.pkl".format(component)))


if __name__ == "__main__":
    main()
