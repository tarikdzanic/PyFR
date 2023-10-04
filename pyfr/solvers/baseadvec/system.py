# -*- coding: utf-8 -*-

from pyfr.solvers.base import BaseSystem
from pyfr.util import memoize


class BaseAdvectionSystem(BaseSystem):
    @memoize
    def _rhs_graphs(self, uinbank, foutbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, foutbank)

        def deps(dk, *names): return self._kdeps(k, dk, *names)

        optimize_memory = self.cfg.getbool('solver', 'optimize-memory', True)

        g1 = self.backend.graph()
        g1.add_mpi_reqs(m['scal_fpts_recv'])

        # Apply positivity-preserving limiter
        g1.add_all(k['eles/limiter'])

        # Compute and store macroscopic state
        g1.add_all(k['eles/macrostate'], deps=k['eles/limiter'])

        # Interpolate the solution to the flux points
        g1.add_all(k['eles/disu'], deps=k['eles/limiter'])

        # Pack and send these interpolated solutions to our neighbours
        g1.add_all(k['mpiint/scal_fpts_pack'], deps=k['eles/disu'])
        for send, pack in zip(m['scal_fpts_send'], k['mpiint/scal_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # Compute the common normal flux at our internal/boundary interfaces
        g1.add_all(k['iint/comm_flux'],
                   deps=k['eles/disu'] + k['mpiint/scal_fpts_pack'])
        g1.add_all(k['bcint/comm_flux'], deps=k['eles/disu'])

        # Make a copy of the solution (if used by source terms)
        g1.add_all(k['eles/copy_soln'], deps=k['eles/limiter'])

        # Interpolate the solution to the quadrature points
        g1.add_all(k['eles/qptsu'], deps=k['eles/limiter'])

        # If separating flux calculation by dimension
        if optimize_memory:
            for i in range(self.ndims):
                if i == 0:
                    tdeps = k['eles/qptsu'] + k['eles/limiter'] + k['eles/copy_soln']
                else:
                    tdeps = k[f'eles/tdivtpcorf_{i-1}']

                # Compute the transformed flux
                g1.add_all(k[f'eles/tdisf_curved_{i}'] + k[f'eles/tdisf_linear_{i}'], deps=tdeps)

                # Compute the transformed divergence of the partially corrected flux
                g1.add_all(k[f'eles/tdivtpcorf_{i}'],
                           deps=k['eles/copy_soln'] + k['eles/disu'] + k['mpiint/scal_fpts_pack'] +
                                k[f'eles/tdisf_curved_{i}'] + k[f'eles/tdisf_linear_{i}'])
        else:
            # Compute the transformed flux
            for l in k['eles/tdisf_curved'] + k['eles/tdisf_linear']:
                g1.add(l, deps=deps(l, 'eles/qptsu', 'eles/limiter'))

            # Compute the transformed divergence of the partially corrected flux
            for l in k['eles/tdivtpcorf']:
                ldeps = deps(l, 'eles/tdisf_curved', 'eles/tdisf_linear',
                            'eles/copy_soln', 'eles/disu')
                g1.add(l, deps=ldeps + k['mpiint/scal_fpts_pack'])
            

        g1.commit()

        g2 = self.backend.graph()

        if optimize_memory:
            g2.add_all(k['eles/copy_soln'])

        # Compute the common normal flux at our MPI interfaces
        g2.add_all(k['mpiint/scal_fpts_unpack'])
        for l in k['mpiint/comm_flux']:
            g2.add(l, deps=deps(l, 'mpiint/scal_fpts_unpack'))

        # Compute the transformed divergence of the corrected flux
        g2.add_all(k['eles/tdivtconf'], deps=k['mpiint/comm_flux'])

        # Obtain the physical divergence of the corrected flux
        for l in k['eles/negdivconf']:
            if optimize_memory:
                g2.add(l, deps=deps(l, 'eles/tdivtconf', 'eles/copy_soln'))
            else:
                g2.add(l, deps=deps(l, 'eles/tdivtconf'))
        g2.commit()

        return g1, g2
