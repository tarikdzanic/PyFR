import numpy as np

from pyfr.solvers.baseadvec import BaseAdvectionElements

class ScalarElements(BaseAdvectionElements):
    privarmap = {2: ['u'],
                 3: ['u']}

    convarmap = {2: ['u'],
                 3: ['u']}

    dualcoeffs = convarmap

    visvarmap = {
        2: [('u', ['u'])],
        3: [('u', ['u'])],
    }

    @staticmethod
    def pri_to_con(pris, cfg):
        return pris

    @staticmethod
    def con_to_pri(cons, cfg):
        return cons

    @staticmethod
    def validate_formulation(ctrl):
        pass

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.tfluxlin')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.limiter')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.elementbounds')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.computeboundselem')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.computeboundsface')

        # Get system parameters
        system = self.cfg.get('solver', 'system')
        if system in ['advection', 'burgers']:
            v = self.cfg.getliteral('solver', 'v')
            assert len(v) == self.ndims
        elif system == 'kpp':
            assert self.ndims == 2
            v = None
        else:
            raise ValueError(f'Unknown system: {system}')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'system': system,
            'v': v
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )
        
        if self.cfg.getbool('solver', 'cbp'):
            tplargs['invvdm'] = self.moninvvdm 
            tplargs['faceinvvdm'] = self.facemoninvvdm 
            tplargs['meanwts'] = self.meanwts
            tplargs['nupts'] = self.nupts
            tplargs['nfpts'] = self.nfpts
            tplargs['mdegs'] = self.basis.ubasis.degrees
            tplargs['nfaces'] = self.nfaces
            tplargs['nfptsperface'] = self.nfptsperface
            tplargs['niters'] = self.cfg.getint('solver', 'niters', 3)

            face_bounds = self.cfg.getbool('solver', 'face-bounds', False)
            elem_bounds = self.cfg.getbool('solver', 'elem-bounds', False)
            glob_bounds = self.cfg.getbool('solver', 'glob-bounds', False)
            assert face_bounds + elem_bounds + glob_bounds == 1, 'Only one bounding method must be enabled.'

            self.kernels['element_bounds'] = lambda uin: self._be.kernel(
                'elementbounds', tplargs=tplargs,
                dims=[self.neles], u=self.scal_upts[uin], x=self.upts_mat,
                bounds=self.bounds, bounds_l=self.bounds_l_int,
                bounds_h=self.bounds_h_int
            )

            if face_bounds:
                self.kernels['compute_bounds'] = lambda : self._be.kernel(
                    'computeboundsface', tplargs=tplargs,
                    dims=[self.neles], uf=self._scal_fpts, xf=self.fpts_mat,
                    bounds=self.bounds
                )
            elif elem_bounds:
                self.kernels['compute_bounds'] = lambda : self._be.kernel(
                    'computeboundselem', tplargs=tplargs,
                    dims=[self.neles], uf=self._scal_fpts, xf=self.fpts_mat,
                    bounds=self.bounds, bounds_l=self.bounds_l_int,
                    bounds_h=self.bounds_h_int
                )
            elif glob_bounds:
                tplargs['global_bounds'] = glob_bounds
                tplargs['gbnds'] = self.cfg.getliteral('solver', 'global-bounds')
                            
            self.kernels['limiter'] = lambda uin: self._be.kernel(
                'limiter', tplargs=tplargs,
                dims=[self.neles], u=self.scal_upts[uin], x=self.upts_mat,
                bounds=self.bounds
            )
