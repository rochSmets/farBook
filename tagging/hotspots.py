#!/usr/bin/env python3

# this is a 2d verison of N hotspots on a plan with elliptical geometry

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics,FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')





numofcells = (200, 200)
mesh = (0.2, 0.2)

Lx = numofcells[0]*mesh[0]
Ly = numofcells[1]*mesh[1]

num_of_spots = 2
fi  = np.radians([20, -20])
psi = np.radians([90, 180])
spot_pos = [[0.5*Lx, 0.0*Ly], [0.5*Lx, 1.0*Ly]]
spot_axis = [[0.2*Lx, 0.2*Ly], [0.2*Lx, 0.2*Ly]]



def polynom(x):
    X = np.fabs(x)
    w = -6*X**5+15*x**4-10*X**3+1
    m = np.clip(np.sign(x), 0, None)

    return m*w


def rotate_coords(pos, beam_id):

    # rotation with fi angle
    wx = (pos[0]-spot_pos[beam_id][0])*np.cos(fi[beam_id])\
        +(pos[1]-spot_pos[beam_id][1])*np.sin(fi[beam_id])
    wy =-(pos[0]-spot_pos[beam_id][0])*np.sin(fi[beam_id])\
        +(pos[1]-spot_pos[beam_id][1])*np.cos(fi[beam_id])

    #rotation with psi angle
    tx = wx
    ty = wy*np.cos(psi[beam_id])

    return [tx, ty]


def density(x, y):

    n = 0.0

    for isp in range(num_of_spots):
        tx, ty = rotate_coords([x, y], isp)
        ux, uy = [tx/spot_axis[isp][0], ty/spot_axis[isp][1]]
        wa = np.sqrt(ux**2+uy**2)

        n += 1.0*polynom(wa)
    #print("x = {0}, y = {1}, n = {2}".format(x, y, n))
    return n


def bx(x, y):
    return 1.0


def by(x, y):
    return 0.0


def bz(x, y):
    return 0.


def v0(x, y):
    return 0.


def vth(x, y):
    return 0.2








def config():

    Simulation(
        smallest_patch_size = 10 ,
        largest_patch_size = 50,
        time_step_nbr = 100,
        final_time = 0.1,
        boundary_types = ["periodic", "periodic"],
        cells=numofcells,
        dl = mesh,
        refinement_boxes = {"L0": {"B0": [(50, 50), (150, 150)]}},
        hyper_resistivity = 0.002,
        resistivity = 0.001,
        diag_options = {"format": "phareh5",
                        "options": {"dir": ".",
                                  "mode":"overwrite"}}
    )



    vMain = {
        "vbulkx": v0, "vbulky": v0, "vbulkz": v0,
        "vthx": vth, "vthy": vth, "vthz": vth,
        "nbr_part_per_cell":100
    }


    MaxwellianFluidModel(
         bx=bx, by=by, bz=bz,
         main={"charge": 1, "density": density, **vMain}
    )


    ElectronModel(closure="isothermal", Te=0.0)


    timestamps = 0.01 * np.arange(10)


    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )

   #for popname in ("protons",):
   #    for name in ["domain", "levelGhost", "patchGhost"]:
   #        ParticleDiagnostics(quantity=name,
   #                            compute_timestamps=timestamps,
   #                            write_timestamps=timestamps,
   #                            population_name=popname)


def main():

    config()

    x = np.arange(20)*0.2
    y = np.arange(20)*0.2
    xv, yv = np.meshgrid(x, y)

    n0 = density(xv, yv)

    # simulator = Simulator(gv.sim)
    # simulator.initialize()
    # simulator.run()


if __name__=="__main__":
    main()
