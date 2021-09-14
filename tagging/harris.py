#!/usr/bin/env python3

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



def config():

    L=0.5

    Simulation(
        smallest_patch_size=10 ,
        largest_patch_size=50,
        time_step_nbr= 100,
        final_time= 0.1,
        boundary_types=["periodic", "periodic"],
        cells=(400,100),
        dl=(0.2, 0.2),
        refinement_boxes={"L0": {"B0": [(150, 2), (250, 48)]}},
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": ".",
                                  "mode":"overwrite"}}
    )

    def density(x, y):
        from pyphare.pharein.global_vars import sim
        Ly = sim.simulation_domain()[1]
        return 0.4 + 1./np.cosh((y-Ly*0.25)/L)**2 + 1./np.cosh((y-Ly*0.75)/L)**2


    def S(y, y0, l):
        return 0.5*(1. + np.tanh((y-y0)/l))


    def by(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.
        dB = 0.1

        x0 = (x - 0.50 * Lx)
        y1 = (y - 0.25 * Ly)
        y2 = (y - 0.75 * Ly)

        dBy1 =  2*dB*x0 * np.exp(-(x0**2 + y1**2)/(sigma)**2)
        dBy2 = -2*dB*x0 * np.exp(-(x0**2 + y2**2)/(sigma)**2)

        return dBy1 + dBy2


    def bx(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.
        dB = 0.1

        x0 = (x - 0.50 * Lx)
        y1 = (y - 0.25 * Ly)
        y2 = (y - 0.75 * Ly)

        dBx1 = -2*dB*y1 * np.exp(-(x0**2 + y1**2)/(sigma)**2)
        dBx2 =  2*dB*y2 * np.exp(-(x0**2 + y2**2)/(sigma)**2)

        v1=-1.
        v2= 1.
        return v1 + (v2-v1)*(S(y,Ly*0.25, 0.5) -S(y, Ly*0.75, 0.5)) + dBx1 + dBx2


    def bz(x, y):
        return 0.


    def b2(x, y):
        return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2


    def T(x, y):
        K = 0.7
        temp = 1./density(x, y)*(K - b2(x, y)*0.5)
        assert np.all(temp >0)
        return temp

    def vx(x, y):
        return 0.


    def vy(x, y):
        return 0.


    def vz(x, y):
        return 0.


    def vthx(x, y):
        return np.sqrt(T(x, y))


    def vthy(x, y):
        return np.sqrt(T(x, y))


    def vthz(x, y):
        return np.sqrt(T(x, y))


    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz,
        "nbr_part_per_cell":100
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density,  **vvv}
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

    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()


if __name__=="__main__":
    main()
