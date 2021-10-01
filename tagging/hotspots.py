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





numofcells = [200, 200]
mesh = [0.2, 0.2]

L = [i*d for i, d in zip(numofcells, mesh)]

num_of_spots = 2
fi  = np.radians([0, 0])
psi = np.radians([40, 40])
spot_pos = [[0.5*L[0], 0.0*L[1]], [0.5*L[0], 1.0*L[1]]]
spot_axis = [[0.4*L[0], 0.4*L[1]], [0.4*L[0], 0.4*L[1]]]



def rect(x):
    return np.where(abs(x)<=1, 1, 0)


def polynom(x):
    X = np.fabs(x)
    w = -6*X**5+15*x**4-10*X**3+1
    return rect(x)*w


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

    x = (0.5+np.arange(numofcells[0]))*mesh[0]
    y = (0.5+np.arange(numofcells[1]))*mesh[1]
    xv, yv = np.meshgrid(x, y)

    n0 = density(xv, yv)

    X = np.arange(numofcells[0]+1)*mesh[0]
    Y = np.arange(numofcells[1]+1)*mesh[1]

    import matplotlib.pyplot as plt
    from matplotlib import rc
    import matplotlib.ticker as ticker

    rc('text', usetex = True)
    rc('font', size=12)
    rc('axes', labelsize='larger')
    rc('mathtext', default='regular')



    fig, ax = plt.subplots(figsize=(6, 5))

    pcm = ax.pcolormesh(X, Y, n0, cmap='viridis_r', edgecolors='face', vmin=0, vmax=1)
    ic = ax.contour(x, y, n0, 8, colors=('k'))

    ax.xaxis.set_major_locator(ticker.LinearLocator(3))
    ax.yaxis.set_major_locator(ticker.LinearLocator(3))

    ax.set_xlabel('$x / l_p$')
    ax.set_ylabel('$y / l_p$')

    plt.title('$\mathrm{Electron \ Density}$')

    cbar = fig.colorbar(pcm, ticks = [0, 0.5, 1], pad = 0.03, aspect = 40)

    plt.savefig('zob.pdf')

    # simulator = Simulator(gv.sim)
    # simulator.initialize()
    # simulator.run()


if __name__=="__main__":
    main()
