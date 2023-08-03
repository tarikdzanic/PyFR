# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)

import numpy as np
from math import gamma as gamma_func

def reflect2D(side, Nx, Ny):
    N = Nx*Ny
    pairs = np.zeros(N, dtype=int)

    for j in range(Ny):
        for i in range(Nx):
            idx = j*Nx + i
            
            if side == 'x':
                ii = Nx - 1 - i
                pidx = j*Nx + ii
            elif side == 'y':
                jj = Ny - 1 - j
                pidx = jj*Nx + i
            
            pairs[idx] = pidx

    return pairs

def reflect3D(side, Nx, Ny, Nz):
    N = Nx*Ny*Nz
    pairs = np.zeros(N, dtype=int)

    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                idx = k*Nx*Ny + j*Nx + i
                
                if side == 'x':
                    ii = Nx - 1 - i
                    pidx = k*Nx*Ny + j*Nx + ii
                elif side == 'y':
                    jj = Ny - 1 - j
                    pidx = k*Nx*Ny + jj*Nx + i
                elif side == 'z':
                    kk = Nz - 1 - k
                    pidx = kk*Nx*Ny + j*Nx + i
                
                pairs[idx] = pidx

    return pairs


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        Ns = [self.cfg.getint('solver', N) for N in ['Nx', 'Ny', 'Nz'][:self.ndims]]
        offsets = np.array([self.cfg.getfloat('solver', off) for off in ['u0', 'v0', 'w0'][:self.ndims]])
        vmax = self.cfg.getfloat('solver', 'vmax')

        # Create velocity bounds
        ubounds = np.zeros((self.ndims, 2))
        ubounds[:, 0] = offsets - vmax
        ubounds[:, 1] = offsets + vmax

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars, rsolver=rsolver,
                             c=self.c, N=Ns, ubounds=ubounds)

class BGKIntInters(TplargsMixin, BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.intcflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )


class BGKMPIInters(TplargsMixin, BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.mpicflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            fl=self._scal_lhs, fr=self._scal_rhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs
        )


class BGKBaseBCInters(TplargsMixin, BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.bgk.kernels.bccflux')

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self.niters = self.cfg.getint('solver', 'niters')
        delta = self.cfg.getint('solver', 'delta')
        Pr = self.cfg.getfloat('solver', 'Pr', 1.0)
        lam = 1.0/gamma_func(delta/2.0) if delta else 0.0

        # Get reflections for wall BCs
        Nx = self.cfg.getint('solver', 'Nx')
        Ny = self.cfg.getint('solver', 'Ny')

        if self.ndims == 2:
            self.Xidxs = reflect2D('x', Nx, Ny)
            self.Yidxs = reflect2D('y', Nx, Ny)
            self.Zidxs = None
        elif self.ndims == 3:
            Nz = self.cfg.getint('solver', 'Nz')
            self.Xidxs = reflect3D('x', Nx, Ny, Nz)
            self.Yidxs = reflect3D('y', Nx, Ny, Nz)
            self.Zidxs = reflect3D('z', Nx, Ny, Nz)

        self._tplargs |= tpl dict(bctype=self.type, niters=self.niters,
                                  pi=np.pi, delta=delta,lam=lam, Pr=Pr,
                                  Xidxs=self.Xidxs, Yidxs=self.Xidxs, Zidxs=self.Xidxs)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, fl=self._scal_lhs,
            magnl=self._mag_pnorm_lhs, nl=self._norm_pnorm_lhs,
            **self._external_vals
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
