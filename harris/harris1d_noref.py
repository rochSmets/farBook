import pyphare.pharein as ph
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
    Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,
        time_step_nbr=100,
        final_time=0.1,
        cells=60,
        dl=0.2,
        resistivity=0.01,
        diag_options={"format": "phareh5",
                      "options": {"dir": "./harris1d_noref",
                                  "mode":"overwrite"}}
    )

    def density(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()[0]
        return 0.5+1./np.cosh((x-L*0.25)/0.5)**2 + 1./np.cosh((x-L*0.75)/0.5)**2

    def S(x, x0, l):
        return 0.5*(1. + np.tanh((x-x0)/l))

    def by(y):
        return  0.;

    def bx(x):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        v1=-1
        v2=1.
        return v1+(v2-v1)*(S(x,Lx*0.25,1)-S(x, Lx*0.75, 1))

    def bz(x):
        return 0.

    def b2(x):
        return bx(x)**2 + by(x)**2 + bz(x)**2

    def T(x):
        K = 1
        return 1./density(x)*(K - b2(x)*0.5)

    def vx(x):
        return 0.

    def vy(x):
        return 0.

    def vz(x):
        return 0.

    def vthx(x):
        return T(x)

    def vthy(x):
        return T(x)

    def vthz(x):
        return T(x)

    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv}
    )

    ElectronModel(closure="isothermal", Te=0.1)

    sim = ph.global_vars.sim

    timestamps = np.arange(0, 0.11, 0.01)

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

def main():
    config()
    simulator = Simulator(gv.sim).initialize().run()
if __name__=="__main__":
    main()

