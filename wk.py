#!/bin/env python3

from dispersion import *
import sys
import os
import numpy as np

def main():
    path = sys.argv[1]
    dx   = float(sys.argv[2])
    dt   = float(sys.argv[3])
    vmin = float(sys.argv[4])
    vmax = float(sys.argv[5])
    component = sys.argv[6]

    B = np.load(os.path.join(path,"{}t.pkl".format(component)))
    Bwk,w, k = fourier(periodize(B), dx, dt)

    print(np.abs(Bwk).min(),np.abs(Bwk).max())

    plot_wk(w, k, Bwk,
            gf,
            vmin=vmin,
            vmax = vmax,
            disp=True,
            kmax = 1,
            wmax = 1.5,
            cmap="Greys",
            color="r",
            filter={"sigma":(0,0)},
            filename=os.path.join(path,"wk.png"))




if __name__ == "__main__":
    main()
