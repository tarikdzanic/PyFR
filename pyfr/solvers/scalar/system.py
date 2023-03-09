from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.scalar.elements import ScalarElements
from pyfr.solvers.scalar.inters import (ScalarIntInters, ScalarMPIInters,
                                       ScalarBaseBCInters)


class ScalarSystem(BaseAdvectionSystem):
    name = 'scalar'

    elementscls = ScalarElements
    intinterscls = ScalarIntInters
    mpiinterscls = ScalarMPIInters
    bbcinterscls = ScalarBaseBCInters
