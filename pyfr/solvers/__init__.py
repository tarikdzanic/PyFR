from pyfr.integrators import get_integrator
from pyfr.solvers.aceuler import ACEulerSystem
from pyfr.solvers.acnavstokes import ACNavierStokesSystem
from pyfr.solvers.base import BaseSystem
from pyfr.solvers.euler import EulerSystem
from pyfr.solvers.navstokes import NavierStokesSystem
from pyfr.solvers.scalar import ScalarSystem
from pyfr.util import subclass_where


def get_solver(backend, rallocs, mesh, initsoln, cfg):
    system = cfg.get('solver', 'system')
    if system in ['advection', 'burgers', 'kpp']:
        systemcls = subclass_where(BaseSystem, name='scalar')
    else:
        systemcls = subclass_where(BaseSystem, name=system)

    # Combine with an integrator to yield the solver
    return get_integrator(backend, systemcls, rallocs, mesh, initsoln, cfg)
