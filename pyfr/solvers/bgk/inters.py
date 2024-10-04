# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)

import functools
import numpy as np

# Generate reflection vectors for specular wall BCs
def reflect(cfg, ndims):
    # Get number of velocity points per dimension
    Ns = [cfg.getint('solver', N) for N in ['Nx', 'Ny', 'Nz'][:ndims]]
    N = np.prod(Ns)

    # Get velocity space indices in grid form
    idxs = np.linspace(0, N-1, N, dtype=int)
    idxsg = np.reshape(idxs, (Ns))

    # Flip velocity indices along normal directions
    out = []
    for i in range(ndims):
        fidxs = np.flip(idxsg, axis=i)
        out.append(np.reshape(fidxs, (-1)))

    return out

class BGKIntInters(BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.intcflux')


        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        delta = self.cfg.getfloat('solver', 'delta', 0.0)
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nuvars=self.nuvars,
                       rsolver=rsolver, c=self.c, delta=delta)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            u=self.umat
        )


class BGKMPIInters(BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.mpicflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        delta = self.cfg.getfloat('solver', 'delta', 0.0)
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nuvars=self.nuvars,
                       rsolver=rsolver, c=self.c, delta=delta)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            u=self.umat
        )


class BGKBaseBCInters(BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.bccflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self.niters = self.cfg.getint('solver', 'niters')
        delta = self.cfg.getfloat('solver', 'delta', 0.0)
        Pr = self.cfg.getfloat('solver', 'Pr', 1.0)

        # Get reflections for wall BCs
        reflidxs = reflect(self.cfg, self.ndims)

        # Get linear system size for DVM
        N = self.ndims + 2

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, nuvars=self.nuvars,
                       c=self.c, u=self.u, bctype=self.type, niters=self.niters,
                       rsolver=rsolver, pi=np.pi, delta=delta, Pr=Pr,
                       reflidxs=reflidxs, N=N)
        
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, fl=self._scal_lhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            u=self.umat, M=self.Mmat, **self._external_vals
        )

class BGKFreeBCInters(BGKBaseBCInters):
    type = 'free'
    cflux_state = 'ghost'

class BGKFixedBCInters(BGKBaseBCInters):
    type = 'fixed'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )

class BGKDiffuseBCInters(BGKBaseBCInters):
    type = 'diffuse'
  
    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c['theta'], = self._eval_opts(['theta'])
        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})

class BGKAdiabaticBCInters(BGKBaseBCInters):
    type = 'adiabatic'
  
    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})

class BGKSpecularBCInters(BGKBaseBCInters):
    type = 'specular'

class BGKInletBCInters(BGKBaseBCInters):
    type = 'inlet'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(
            ['rho', 'u', 'v', 'w'][:self.ndims + 1], lhs
        )

class BGKOutletBCInters(BGKBaseBCInters):
    type = 'outlet'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['p'], lhs)

class BGKPressureThetaBCInters(BGKBaseBCInters):
    type = 'prth'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['p'], lhs)
        self.c |= self._exp_opts(['theta'], lhs)
