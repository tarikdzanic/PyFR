from pyfr.solvers.base import BaseSystem
from pyfr.util import memoize


class BaseAdvectionSystem(BaseSystem):
    @memoize
    def _rhs_graphs(self, uinbank, foutbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, foutbank)

        def deps(dk, *names): return self._kdeps(k, dk, *names)

        g1 = self.backend.graph()
        g1.add_mpi_reqs(m['scal_fpts_recv'])

        # # Interpolate the solution to the flux points
        g1.add_all(k['eles/disu'])

        # Pack and send these interpolated solutions to our neighbours
        g1.add_all(k['mpiint/scal_fpts_pack'], deps=k['eles/disu'])
        for send, pack in zip(m['scal_fpts_send'], k['mpiint/scal_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # Compute the common normal flux at our internal/boundary interfaces
        g1.add_all(k['iint/comm_flux'],
                   deps=k['eles/disu'] + k['mpiint/scal_fpts_pack'])
        g1.add_all(k['bcint/comm_flux'], deps=k['eles/disu'])

        # Make a copy of the solution (if used by source terms)
        g1.add_all(k['eles/copy_soln'])

        # Interpolate the solution to the quadrature points
        g1.add_all(k['eles/qptsu'])

        # Compute the transformed flux
        for l in k['eles/tdisf_curved'] + k['eles/tdisf_linear']:
            ldeps = deps(l, 'eles/qptsu')
            g1.add(l, deps=ldeps)

        # Compute the transformed divergence of the partially corrected flux
        for l in k['eles/tdivtpcorf']:
            ldeps = deps(l, 'eles/tdisf_curved', 'eles/tdisf_linear',
                         'eles/copy_soln', 'eles/disu')
            g1.add(l, deps=ldeps + k['mpiint/scal_fpts_pack'])
        g1.commit()

        g2 = self.backend.graph()

        # Compute the common normal flux at our MPI interfaces
        g2.add_all(k['mpiint/scal_fpts_unpack'])
        for l in k['mpiint/comm_flux']:
            g2.add(l, deps=deps(l, 'mpiint/scal_fpts_unpack'))

        # Compute the transformed divergence of the corrected flux
        g2.add_all(k['eles/tdivtconf'], deps=k['mpiint/comm_flux'])

        # Obtain the physical divergence of the corrected flux
        for l in k['eles/negdivconf']:
            g2.add(l, deps=deps(l, 'eles/tdivtconf'))
        g2.commit()

        return g1, g2

    @memoize
    def _preproc_graphs(self, uinbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, None)

        def deps(dk, *names): return self._kdeps(k, dk, *names)
        bound_method = self.cfg.get('solver', 'bounds', None)

        if bound_method == 'local-face':
            raise NotImplementedError()
            g1 = self.backend.graph()

            # Interpolate the solution to the flux points
            if 'eles/element_bounds' in k:
                g1.add_all(k['eles/disu'])

            # Compute local bounds within element
            g1.add_all(k['eles/element_bounds'])

            # Pack and send the bounds values to neighbors
            g1.add_all(k['mpiint/bounds_fpts_pack'], deps=k['eles/disu'])
            for send, pack in zip(m['bounds_fpts_send'], k['mpiint/bounds_fpts_pack']):
                g1.add_mpi_req(send, deps=[pack])

            # Compute common bounds at internal/boundary interfaces
            g1.add_all(k['iint/comm_bounds'], deps=k['eles/disu'])
            g1.add_all(k['bcint/comm_bounds'], deps=k['eles/disu'])

            if 'mpiint/comm_bounds' in k:
                # Compute common entropy minima at MPI interfaces
                g2 = self.backend.graph()

                g2.add_all(k['mpiint/bounds_fpts_unpack'])
                for l in k['mpiint/comm_bounds']:
                    g2.add(l, deps=deps(l, 'mpiint/bounds_fpts_unpack'))

                g2.add_all(k['eles/compute_bounds'], deps=k['mpiint/comm_bounds'])
                g2.commit()

                return g1, g2
            else:
                g1.add_all(k['eles/compute_bounds'], deps=k['eles/element_bounds'] +
                                                          k['iint/comm_bounds'] +
                                                          k['bcint/comm_bounds'])
                g1.commit()
                return g1,
        if bound_method == 'local-element':
            g1 = self.backend.graph()
            g1.add_mpi_reqs(m['bounds_l_fpts_recv'])
            g1.add_mpi_reqs(m['bounds_h_fpts_recv'])
            g1.add_mpi_reqs(m['bounds_e_fpts_recv'])

            # Interpolate the solution to the flux points
            if 'eles/element_bounds' in k:
                g1.add_all(k['eles/disu'])

            # Compute local bounds within element
            g1.add_all(k['eles/element_bounds'])

            # Pack and send the bounds values to neighbors
            g1.add_all(k['mpiint/bounds_l_fpts_pack'], deps=k['eles/element_bounds'])
            for send, pack in zip(m['bounds_l_fpts_send'], k['mpiint/bounds_l_fpts_pack']):
                g1.add_mpi_req(send, deps=[pack])
            g1.add_all(k['mpiint/bounds_h_fpts_pack'], deps=k['eles/element_bounds'])
            for send, pack in zip(m['bounds_h_fpts_send'], k['mpiint/bounds_h_fpts_pack']):
                g1.add_mpi_req(send, deps=[pack])
            g1.add_all(k['mpiint/bounds_e_fpts_pack'], deps=k['eles/element_bounds'])
            for send, pack in zip(m['bounds_e_fpts_send'], k['mpiint/bounds_e_fpts_pack']):
                g1.add_mpi_req(send, deps=[pack])

            # Compute common entropy minima at internal/boundary interfaces
            g1.add_all(k['iint/comm_bounds_l'], deps=k['eles/element_bounds'])
            g1.add_all(k['iint/comm_bounds_h'], deps=k['eles/element_bounds'])
            g1.add_all(k['iint/comm_bounds_e'], deps=k['eles/element_bounds'])
            g1.add_all(k['bcint/comm_bounds_l'],
                    deps=k['eles/element_bounds'] + k['eles/disu'])
            g1.add_all(k['bcint/comm_bounds_h'],
                    deps=k['eles/element_bounds'] + k['eles/disu'])
            g1.add_all(k['bcint/comm_bounds_e'],
                    deps=k['eles/element_bounds'] + k['eles/disu'])

            if 'mpiint/comm_bounds_l' in k:
                # Compute common entropy minima at MPI interfaces
                g2 = self.backend.graph()

                g2.add_all(k['mpiint/bounds_l_fpts_unpack'])
                for l in k['mpiint/comm_bounds_l']:
                    g2.add(l, deps=deps(l, 'mpiint/bounds_l_fpts_unpack'))
                g2.add_all(k['mpiint/bounds_h_fpts_unpack'])
                for l in k['mpiint/comm_bounds_h']:
                    g2.add(l, deps=deps(l, 'mpiint/bounds_h_fpts_unpack'))
                g2.add_all(k['mpiint/bounds_e_fpts_unpack'])
                for l in k['mpiint/comm_bounds_e']:
                    g2.add(l, deps=deps(l, 'mpiint/bounds_e_fpts_unpack'))

                g2.add_all(k['eles/compute_bounds'], deps=k['mpiint/comm_bounds_l'] +
                                                          k['mpiint/comm_bounds_h'] +
                                                          k['mpiint/comm_bounds_e'])
                g2.commit()

                return g1, g2
            else:
                g1.add_all(k['eles/compute_bounds'], deps=k['iint/comm_bounds_l'] +
                                                          k['iint/comm_bounds_h'] +
                                                          k['iint/comm_bounds_e'] +
                                                          k['bcint/comm_bounds_l'] +
                                                          k['bcint/comm_bounds_h'] +
                                                          k['bcint/comm_bounds_e'])
                g1.commit()
                return g1,
        else:
            return []

    def postproc(self, uinbank):
        k, _ = self._get_kernels(uinbank, None)

        self.backend.run_kernels(k['eles/limiter'])
