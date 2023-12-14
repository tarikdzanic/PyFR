from pyfr.solvers.base import BaseElements
from pyfr.quadrules import get_quadrule

import numpy as np

class BaseAdvectionElements(BaseElements):
    @property
    def _scratch_bufs(self):
        if 'flux' in self.antialias:
            bufs = {'scal_fpts', 'scal_qpts', 'vect_qpts'}
        else:
            bufs = {'scal_fpts', 'vect_upts'}

        bufs |= {'scal_upts_cpy'}

        return bufs

    def set_backend(self, backend, nscalupts, nonce, linoff):
        super().set_backend(backend, nscalupts, nonce, linoff)

        kernels = self.kernels

        # Register pointwise kernels with the backend
        self._be.pointwise.register(
            'pyfr.solvers.baseadvec.kernels.negdivconf'
        )

        # What anti-aliasing options we're running with
        fluxaa = 'flux' in self.antialias

        # What the source term expressions (if any) are a function of
        plocsrc = self._ploc_in_src_exprs
        solnsrc = self._soln_in_src_exprs

        # Source term kernel arguments
        srctplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'srcex': self._src_exprs
        }

        # Interpolation from elemental points
        kernels['disu'] = lambda uin: self._be.kernel(
            'mul', self.opmat('M0'), self.scal_upts[uin],
            out=self._scal_fpts
        )

        if fluxaa and self.basis.order > 0:
            kernels['qptsu'] = lambda uin: self._be.kernel(
                'mul', self.opmat('M7'), self.scal_upts[uin],
                out=self._scal_qpts
            )

        # First flux correction kernel
        if fluxaa and self.basis.order > 0:
            kernels['tdivtpcorf'] = lambda fout: self._be.kernel(
                'mul', self.opmat('(M1 - M3*M2)*M9'), self._vect_qpts,
                out=self.scal_upts[fout]
            )
        elif self.basis.order > 0:
            kernels['tdivtpcorf'] = lambda fout: self._be.kernel(
                'mul', self.opmat('M1 - M3*M2'), self._vect_upts,
                out=self.scal_upts[fout]
            )

        # Second flux correction kernel
        kernels['tdivtconf'] = lambda fout: self._be.kernel(
            'mul', self.opmat('M3'), self._scal_fpts,
            out=self.scal_upts[fout], beta=float(self.basis.order > 0)
        )

        # Transformed to physical divergence kernel + source term
        plocupts = self.ploc_at('upts') if plocsrc else None
        solnupts = self._scal_upts_cpy if solnsrc else None

        kernels['copy_soln'] = lambda uin: self._be.kernel(
            'copy', self._scal_upts_cpy, self.scal_upts[uin]
        )
        kernels['copy_soln2'] = lambda uin: self._be.kernel(
            'copy', self.scal_upts[uin], self._scal_upts_cpy
        )

        kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconf', tplargs=srctplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, u=solnupts
        )

        if self.cfg.getbool('solver', 'cbp') and self.basis.order > 0:
            # Create monomial VDM
            degs = self.basis.ubasis.degrees
            upts = self.basis.upts
            V = np.ones((self.nupts, self.nupts))

            for i in range(self.nupts):
                for j, dd in enumerate(degs[i]):
                    V[i,:] *= upts[:,j]**dd
            
            self.moninvvdm = np.linalg.inv(V.T)
            ub = self.basis.ubasis
            self.meanwts = ub.invvdm[:,0]/np.sum(ub.invvdm[:,0])

            # npts = 4*self.basis.order + 1
            # ux = np.linspace(-1, 1, npts)
            # uxx, uyy = np.meshgrid(ux, ux, indexing='xy')
            # ug = np.zeros((npts**2, 2))
            # ug[:,0] = np.reshape(uxx, (-1))
            # ug[:,1] = np.reshape(uyy, (-1))
            self.upts_mat = self._be.const_matrix(upts)
            assert self.ndims == 2, "Newton's method and face search only implemented for 2D."

            # Setup LMP bounds
            self.nfaces = len(self.nfacefpts)
            assert all(nf == self.nfacefpts[0] for nf in self.nfacefpts), 'All faces must have same number of fpts.'
            self.nfptsperface = self.nfacefpts[0]
            rule = self.cfg.get('solver-interfaces-line', 'flux-pts')
            npts = self.basis.npts_for_face['line'](self.basis.order)
            x_fpts = np.array(get_quadrule('line', rule, npts).pts)
            # Create monomial VDM for face solution
            V = np.empty((len(x_fpts), len(x_fpts)))
            for i in range(len(x_fpts)):
                V[i,:] = x_fpts**i
            self.facemoninvvdm = np.linalg.inv(V.T)
            self.fpts_mat = self._be.const_matrix(np.atleast_2d(x_fpts).T)
            self.bounds = self._be.matrix((3, self.neles),
                                           extent=nonce + 'bounds',
                                           tags={'align'})

            self.bounds_l_int = self._be.matrix((self.nfaces, self.neles),
                                                extent=nonce + 'bounds_l_int',
                                                tags={'align'})
            self.bounds_h_int = self._be.matrix((self.nfaces, self.neles),
                                                extent=nonce + 'bounds_h_int',
                                                tags={'align'})
            self.bounds_e_int = self._be.matrix((self.nfaces, self.neles),
                                                extent=nonce + 'bounds_e_int',
                                                tags={'align'})

        # In-place solution filter
        if self.cfg.getint('soln-filter', 'nsteps', '0'):
            def modal_filter(uin):
                mul = self._be.kernel(
                    'mul', self.opmat('M10'), self.scal_upts[uin],
                    out=self._scal_upts_temp
                )
                copy = self._be.kernel(
                    'copy', self.scal_upts[uin], self._scal_upts_temp
                )

                return self._be.ordered_meta_kernel([mul, copy])

            kernels['modal_filter'] = modal_filter

        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')
        if shock_capturing == 'entropy-filter':
            tags = {'align'}

            # Allocate one minimum entropy value per interface
            self.nfaces = len(self.nfacefpts)
            ext = nonce + 'entmin_int'
            self.entmin_int = self._be.matrix((self.nfaces, self.neles),
                                              tags=tags, extent=ext)

            # Setup nodal/modal operator matrices
            self.vdm = self._be.const_matrix(self.basis.ubasis.vdm.T)
            self.invvdm = self._be.const_matrix(self.basis.ubasis.invvdm.T)
        else:
            self.entmin_int = None

    def get_entmin_int_fpts_for_inter(self, eidx, fidx):
        return (self.entmin_int.mid,), (fidx,), (eidx,)

    def get_entmin_bc_fpts_for_inter(self, eidx, fidx):
        nfp = self.nfacefpts[fidx]
        return (self.entmin_int.mid,)*nfp, (fidx,)*nfp, (eidx,)*nfp

    def get_bounds_l_int_fpts_for_inter(self, eidx, fidx):
        return (self.bounds_l_int.mid,), (fidx,), (eidx,)

    def get_bounds_l_bc_fpts_for_inter(self, eidx, fidx):
        nfp = self.nfacefpts[fidx]
        return (self.bounds_l_int.mid,)*nfp, (fidx,)*nfp, (eidx,)*nfp

    def get_bounds_h_int_fpts_for_inter(self, eidx, fidx):
        return (self.bounds_h_int.mid,), (fidx,), (eidx,)

    def get_bounds_h_bc_fpts_for_inter(self, eidx, fidx):
        nfp = self.nfacefpts[fidx]
        return (self.bounds_h_int.mid,)*nfp, (fidx,)*nfp, (eidx,)*nfp

    def get_bounds_e_int_fpts_for_inter(self, eidx, fidx):
        return (self.bounds_e_int.mid,), (fidx,), (eidx,)

    def get_bounds_e_bc_fpts_for_inter(self, eidx, fidx):
        nfp = self.nfacefpts[fidx]
        return (self.bounds_e_int.mid,)*nfp, (fidx,)*nfp, (eidx,)*nfp
