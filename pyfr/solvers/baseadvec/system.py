# -*- coding: utf-8 -*-

from pyfr.solvers.base import BaseSystem


class BaseAdvectionSystem(BaseSystem):
    _nqueues = 2

    def rhs(self, t, uinbank, foutbank):
        self.rhs_LO(t, uinbank, foutbank)

    def rhs_LO(self, t, uinbank, foutbank):
        runall = self.backend.runall
        q1, q2 = self._queues
        kernels = self._kernels

        self._bc_inters.prepare(t)

        self.eles_scal_upts_inb.active = uinbank
        self.eles_scal_upts_outb.active = foutbank

        self.calc_LOupwind_divf(runall, q1, q2, kernels, t)
        # self.calc_LOcentered_divf(runall, q1, q2, kernels, t)
        # self.calc_HO_divf(runall, q1, q2, kernels, t)
        self.calc_LO_residual(runall, q1, q2, kernels, t)

        # self.calc_HO_divf(runall, q1, q2, kernels, t)
        # self.calc_HO_residual(runall, q1, q2, kernels, t)

        self.calc_RD_divf(runall, q1, q2, kernels, t)

    def calc_LOupwind_divf(self, runall, q1, q2, kernels, t):
        # Calculate low-order gradients
        q1 << kernels['eles', 'disu_LO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_LO_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()

        # q1 << kernels['eles', 'tdisf']()
        q1 << kernels['eles', 'divf_LO']()
        # q1 << kernels['eles', 'tdivtpcorf_LOsubcell']()

        q1 << kernels['iint', 'comm_flux_centered']()
        q1 << kernels['bcint', 'comm_flux_centered'](t=t)

        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux_centered']()

        q1 << kernels['eles', 'tdivtconf_LO']()


    def calc_LOcentered_divf(self, runall, q1, q2, kernels, t):
        # Calculate low-order gradients
        q1 << kernels['eles', 'disu_LO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_LO_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()

        q1 << kernels['eles', 'tdisf']()
        q1 << kernels['eles', 'tdivtpcorf_LOcentered']()

        q1 << kernels['iint', 'comm_flux_centered']()
        q1 << kernels['bcint', 'comm_flux_centered'](t=t)

        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux_centered']()
        q1 << kernels['eles', 'tdivtconf_LO']()

    def calc_HO_divf(self, runall, q1, q2, kernels, t):
        # Calculate low-order gradients
        q1 << kernels['eles', 'disu_HO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_HO_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()

        q1 << kernels['eles', 'tdisf']()
        q1 << kernels['eles', 'tdivtpcorf_HO']()

        q1 << kernels['iint', 'comm_flux']()
        q1 << kernels['bcint', 'comm_flux'](t=t)

        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux']()
        q1 << kernels['eles', 'tdivtconf_HO']()

    def calc_LO_residual(self, runall, q1, q2, kernels, t):
        q1 << kernels['eles', 'residual']()
        # Map u to f shape and calculate gradient
        q1 << kernels['eles', 'disu_LO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_LO_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()

        q1 << kernels['eles', 'tdisu']()
        q1 << kernels['eles', 'tdivtpcorf_LO']()

        q1 << kernels['iint', 'comm_sol_centered']()
        q1 << kernels['bcint', 'comm_sol_centered'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_sol_centered']()
        q1 << kernels['eles', 'tdivtconf_LO']()
        q1 << kernels['eles', 'normalizeresidual']()

    def calc_HO_residual(self, runall, q1, q2, kernels, t):
        q1 << kernels['eles', 'residual']()
        # Map u to f shape and calculate gradient
        q1 << kernels['eles', 'disu_HO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_HO_int']()
        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()
        if ('eles', 'copy_soln_at_fpts') in kernels:
            q1 << kernels['eles', 'copy_soln_at_fpts']()

        q1 << kernels['eles', 'tdisu']()
        q1 << kernels['eles', 'tdivtpcorf_HO']()

        q1 << kernels['iint', 'comm_sol_centered']()
        q1 << kernels['bcint', 'comm_sol_centered'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_sol_centered']()
        q1 << kernels['eles', 'tdivtconf_HO']()
        q1 << kernels['eles', 'normalizeresidual']()

    def calc_RD_divf(self, runall, q1, q2, kernels, t):
        # Calculate LO comm_flux and set it to _scal_fpts_cpy
        q1 << kernels['eles', 'disu_LO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_LO_int']()

        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()

        q1 << kernels['iint', 'comm_flux']()
        q1 << kernels['bcint', 'comm_flux'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])

        q1 << kernels['mpiint', 'comm_flux']()
        q1 << kernels['eles', 'copy_soln_at_fpts']()
        runall([q1])

        # Calculate HO comm_flux and set it to _scal_fpts
        q1 << kernels['eles', 'disu_HO_ext']()
        q1 << kernels['mpiint', 'scal_fpts_pack']()
        runall([q1])
        q1 << kernels['eles', 'disu_HO_int']()
        q1 << kernels['eles', 'limitinterp']()

        if ('eles', 'copy_soln') in kernels:
            q1 << kernels['eles', 'copy_soln']()

        q1 << kernels['iint', 'comm_flux']()
        q1 << kernels['bcint', 'comm_flux'](t=t)
        q2 << kernels['mpiint', 'scal_fpts_send']()
        q2 << kernels['mpiint', 'scal_fpts_recv']()
        q2 << kernels['mpiint', 'scal_fpts_unpack']()
        runall([q1, q2])
        q1 << kernels['mpiint', 'comm_flux']()

        # Calculate RD divf without interface contributions
        q1 << kernels['eles', 'riemanndifference']()

        # Choose flux based on alpha 
        q1 << kernels['eles', 'blendintflux']()
        q1 << kernels['eles', 'tdivtconf_RD']()

        if ('eles', 'tdivf_qpts') in kernels:
            raise NotImplementedError()
        else:
            q1 << kernels['eles', 'negdivconf'](t=t)
        runall([q1])
