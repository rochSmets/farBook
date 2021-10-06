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

# num_of_spots = 2
# fi  = np.radians([0, 0])
# psi = np.radians([-40, 40])
# spot_pos = [[0.5*L[0], 0.0*L[1]], [0.5*L[0], 1.0*L[1]]]
# spot_axis = [[0.1*L[0], 0.5*L[1]], [0.1*L[0], 0.5*L[1]]]
# spot_width = [0.2, 0.2]

num_of_spots = 1
fi  = np.radians([90])
psi = np.radians([40])
spot_pos = [[0.5*L[0], 0.5*L[1]]]
spot_axis = [[0.3*L[0], 0.3*L[1]]]
spot_width = [[0.5]]



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
    ty = wy/np.cos(psi[beam_id])

    return [tx, ty]


def density(x, y):
    n = 0.0

    for isp in range(num_of_spots):
        tx, ty = rotate_coords([x, y], isp)
        ux, uy = [tx/spot_axis[isp][0], ty/spot_axis[isp][1]]
        ut = np.sqrt(ux**2+uy**2)

        n += 1.0*polynom(ut)
    return n


def b_and_u(x, y, isp):
    u = [None]*3

    tx, ty = rotate_coords([x, y], isp)
    ux, uy = [tx/spot_axis[isp][0], ty/spot_axis[isp][1]]
    ut = np.clip(np.sqrt(ux**2+uy**2), 0.001, None)
    wt = (ut-1)/spot_width[isp]

    vx = +uy/spot_axis[isp][1]
    vy = -ux/spot_axis[isp][0]
    vt = np.clip(np.sqrt(vx**2+vy**2), 0.001, None)

    # this normalization term is mandatory to ensure div B = 0 in an ellipsis
    z = spot_axis[isp][1]*vt/ut

    u[0] = (vx/vt)*np.cos(fi[isp])-(vy/vt)*np.sin(fi[isp])*np.cos(psi[isp])
    u[1] = (vx/vt)*np.sin(fi[isp])+(vy/vt)*np.cos(fi[isp])*np.cos(psi[isp])
    u[2] =                         (vy/vt)*               np.sin(psi[isp])

    return 1.0*polynom(wt)*z, u


def bx(x, y):
    bx = 0.0
    for isp in range(num_of_spots):
        b, u = b_and_u(x, y, isp)
        bx += b*u[0]
    return bx


def by(x, y):
    by = 0.0
    for isp in range(num_of_spots):
        b, u = b_and_u(x, y, isp)
        by += b*u[1]
    return by


def bz(x, y):
    bz = 0.0
    for isp in range(num_of_spots):
        b, u = b_and_u(x, y, isp)
        bz += b*u[2]
    return bz


def v0(x, y):
    return 0.


def vth(x, y):
    return np.sqrt(0.2)


#/home/smets/codeS/pYwi/pywi/runs/run.py


def az(x, y):
    from numpy.fft import fft2
    from numpy.fft import ifft2

    N = [n+1 for n in numofcells]

    x = np.linspace(0, +L[0], N[0])
    y = np.linspace(0, +L[1], N[1])

    xv, yv = np.meshgrid(x, y)

    bx_ = bx(xv, yv)+1j*np.zeros(yv.shape)
    by_ = by(xv, yv)+1j*np.zeros(xv.shape)

    def k_(axis, L, N):
        N_ = N[axis]
        L_ = L[axis]
        m = np.linspace(0, N_-1, N_)
        mplus  =  m[            :int(N_/2)+1   ]
        mminus = -m[int((N_-1)/2):0          :-1]
        return np.concatenate((mplus, mminus))*2*np.pi/L_

    kx = k_(0, L, N)
    ky = k_(1, L, N)

    kxv, kyv = np.meshgrid(kx, ky)

    k2 = np.square(kxv)+np.square(kyv)
    k2[0][0] = 1.0

    BX = fft2(bx_)
    BY = fft2(by_)

    AZ = np.divide(1j*(kxv*BY-kyv*BX), k2)
    AZ[0][0] = 0.0

    az = ifft2(AZ).real

    return az[:-1, :-1]



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

    X = np.arange(numofcells[0]+1)*mesh[0]
    Y = np.arange(numofcells[1]+1)*mesh[1]

    import matplotlib.pyplot as plt
    from matplotlib import rc
    import matplotlib.ticker as ticker

    rc('text', usetex = True)
    rc('font', size=12)
    rc('axes', labelsize='larger')
    rc('mathtext', default='regular')

    a = az(xv, yv)
    #   z = density(xv, yv)
    #   z, u = b_and_u(xv, yv, 0)
    z = bz(xv, yv)


    Bx = bx(xv, yv)
    By = by(xv, yv)

    dbxdx = np.gradient(Bx, axis=1)
    dbydy = np.gradient(By, axis=0)
    divB = dbxdx+dbydy
    print("       div B = {0:4.6f}".format(np.fabs(dbxdx+dbydy).max()))
    print("pseudo-div B = {0:4.6f}".format(np.fabs(dbxdx+dbydy).max()*2/np.fabs(dbxdx-dbydy).max()))

    minmax=[-0.4, +0.4]
    fig, ax = plt.subplots(figsize=(6, 5))

    # pcm = ax.pcolormesh(X, Y, z, cmap='viridis_r', edgecolors='face')
    pcm = ax.pcolormesh(X, Y, z, cmap='viridis_r', edgecolors='face', vmin=minmax[0], vmax=minmax[1])
    ic = ax.contour(x, y, a, 12, colors=('k'), linewidths=(0.2))

    ax.xaxis.set_major_locator(ticker.LinearLocator(3))
    ax.yaxis.set_major_locator(ticker.LinearLocator(3))

    ax.set_xlabel('$x / l_p$')
    ax.set_ylabel('$y / l_p$')

    plt.title('$\mathrm{2\ Hot \ Spots}$')

    cbar = fig.colorbar(pcm, ticks = [minmax[0], minmax[1]], pad = 0.03, aspect = 40)

    plt.savefig('zob.png')

    # simulator = Simulator(gv.sim)
    # simulator.initialize()
    # simulator.run()


if __name__=="__main__":
    main()
