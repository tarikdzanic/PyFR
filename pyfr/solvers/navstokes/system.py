from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.navstokes.elements import NavierStokesElements
from pyfr.solvers.navstokes.inters import (NavierStokesBaseBCInters,
                                           NavierStokesIntInters,
                                           NavierStokesMPIInters,
                                           NavierStokesPintInters)


class NavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'navier-stokes'

    elementscls = NavierStokesElements
    intinterscls = NavierStokesIntInters
    mpiinterscls = NavierStokesMPIInters
    pintinterscls = NavierStokesPintInters
    bbcinterscls = NavierStokesBaseBCInters
