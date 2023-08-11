# -*- coding: utf-8 -*-

from ctypes.wintypes import PSIZE
from pyfr.solvers.baseadvec import BaseAdvectionElements

import functools
from math import gamma as gamma_func
import numpy as np

# Setup velocity space and integrator
def setup_BGK(cfg, ndims):
    # Get number of velocity points and velocity offset per dimension
    Ns = [cfg.getint('solver', N) for N in ['Nx', 'Ny', 'Nz'][:ndims]]
    offsets = np.array([cfg.getfloat('solver', off) for off in ['u0', 'v0', 'w0'][:ndims]])
    vmax = cfg.getfloat('solver', 'vmax')

    # Create velocity bounds
    mins = list(offsets - vmax)
    maxs = list(offsets + vmax)

    # Helper function to create 1D trapezoidal rule
    linwts = lambda N, mass: (np.array([0.5] + list(np.ones(N)[1:-1]) + [0.5]))*mass/(N-1)

    # Create velocity/integrator grid
    ug = np.meshgrid(*[np.linspace(ul, uh, N) for ul, uh, N in zip(mins, maxs, Ns)], indexing='ij')
    Mg = functools.reduce(np.multiply, np.ix_(*[linwts(N, uh - ul) for ul, uh, N in zip(mins, maxs, Ns)]))

    # Reduce grids to 1D and 2D arrays
    M = Mg.reshape(-1)
    u = np.empty((np.prod(Ns), len(Ns)))
    for i in range(len(Ns)):
        u[:,i] = ug[i].reshape(-1)

    # Compute collision invariants
    psi = np.empty((len(u), ndims+2))
    psi[:,0] = 1.0
    for i in range(ndims):
        psi[:,i+1] = u[:,i]
    psi[:,-1] = 0.5*np.linalg.norm(u, axis=1)**2

    return [u, M, psi]

# Computes discretely conservative Maxwellian for a macroscopic solution U using discrete velocity model
def iterate_DVM(U, u, ndims, psi, M, gamma, niters, delta):
    def compute_discrete_maxwellian(alpha):
        # Compute the macro/micro velocity defect
        dv2 = 0
        for i in range(ndims):
            dv2 += (u[...,i] - alpha[i+2])**2

        # Compute Maxwellian (monatomic)
        M = (alpha[0]*np.exp(-alpha[1]*dv2))

        return M

    # Change local variables into alpha vector
    rho, E = U[0], U[-1]
    vs = [rhov/rho for rhov in U[1:-1]]
    theta = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in vs))/rho

    alpha = np.empty(ndims+2)
    alpha[0] = rho/(2*np.pi*theta)**(ndims/2.0)
    alpha[1] = 1.0/(2*theta)
    for i in range(ndims):
        alpha[i+2] = vs[i]

    # Compute initial guess for Maxwellian (analytic distribution)
    g = compute_discrete_maxwellian(alpha)

    # Perform Newton iterations to find optimal Maxwellian
    F = [0.0]*(ndims+2)
    for _ in range(niters):
        # Derivatives with respect to alpha
        Q = [None]*(ndims+2)
        Q[0] = 1.0/alpha[0]
        Q[1] = 0.0

        for i in range(ndims):
            Q[1] += -(u[...,i] - alpha[i+2])**2
            Q[i+2] = 2*alpha[1]*(u[...,i] - alpha[i+2])

        # Compute Jacobian
        J = np.zeros((ndims+2, ndims+2))
        for ivar in range(ndims+2):
            psig = psi[:,ivar]*g
            F[ivar] = np.dot(M, psig) - U[ivar]

            for jvar in range(ndims+2):
                J[ivar, jvar] = np.dot(M, Q[jvar]*psig)

        # Take Newton step and compute new discrete Maxwellian
        alpha = alpha - np.linalg.solve(J, F)
        g = compute_discrete_maxwellian(alpha)

    # Print warning if macroscopic state residual exceeds 1e-6
    if np.max(F) > 1e-6:
        print(f'DVM did not converge for solution {U} with residual {F}.')

    return g

class BGKElements(BaseAdvectionElements):    
    formulations = ['std', 'dual']

    privarmap = {2: ['rho', 'u', 'v', 'p'],
                 3: ['rho', 'u', 'v', 'w', 'p']}

    privarmap2 = {2: ['rho', 'u', 'v', 'p', 'sxy'],
                  3: ['rho', 'u', 'v', 'w', 'p', 'sxy', 'sxz', 'syz']}

    convarmap = {2: ['1'],
                 3: ['1']}

    dualcoeffs = convarmap

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('pressure', ['p']),
            ('strain', ['sxy'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('pressure', ['p']),
            ('strain', ['sxy', 'sxz', 'syz'])]
    }


    def __init__(self, basiscls, eles, cfg):
        self.ndims = eles.shape[2]

        [self.u, self.M, self.psi] = setup_BGK(cfg, self.ndims)
        self.delta = cfg.getint('solver', 'delta')
        self.iterate_ICs = cfg.getbool('solver', 'iterate_ICs', True)

        self.nuvars = len(self.u)
        self.nvars = 2*self.nuvars if self.delta else self.nuvars
        self.nmvars = self.ndims + 2

        super().__init__(basiscls, eles, cfg)

    # Compute Maxwellian state from primitive initial conditions
    def pri_to_con(self, pris, cfg):
        # Convert primitive macroscopic state to conserved macroscopic state
        cons = self.macropri_to_macrocon(pris, cfg)

        # Allocate initial distribution function
        f = np.zeros((self.nupts, self.nuvars, self.neles))
        gamma = cfg.getfloat('constants', 'gamma')

        niters = cfg.getint('solver', 'niters') if self.iterate_ICs else 0 # Large iteration count for ICs
        for uidx in range(self.nupts):
            for eidx in range(self.neles):
                # Get local conserved state variables
                cons_local = np.zeros(self.ndims+2)
                for i in range(self.ndims+2):
                    cons_local[i] = cons[i] if np.isscalar(cons[i]) else cons[i][uidx, eidx]

                # Compute local Maxwellian
                f[uidx, :, eidx] = iterate_DVM(cons_local, self.u, self.ndims, self.psi, self.M, gamma, niters, self.delta)

        if self.delta:
            fg = np.zeros((self.nupts, self.nvars, self.neles))
            fg[:,:self.nuvars,:] = f
            return fg
        else:
            return f

    # Compute macroscopic primitive state variables (moments) from distribution function
    @staticmethod
    def con_to_pri(f, cfg, M, u, psi, ndims):
        pris = []
        for i in range(ndims+2):
            pris.append(np.einsum('i,ijk->jk', M*psi[...,i], f))

        return pris
    
    @staticmethod
    def con_to_vis(f, cfg, M, u, psi, ndims):
        # Compute primitive variables
        pris = BGKElements.con_to_pri(f, cfg, M, u, psi, ndims)

        # Compute and append off-diagonal molecular stresses
        if ndims == 2:
            psi2 = [u[:,0]*u[:,1]]
        elif ndims == 3:
            psi2 = [u[:,0]*u[:,1], u[:,0]*u[:,2], u[:,1]*u[:,2]]
 
        for i in range(len(psi2)):
            pris.append(np.einsum('i,ijk->jk', M*psi2[i], f))

        return pris

    @staticmethod
    def macrocon_to_macropri(cons, cfg):
        rho, E = cons[0], cons[-1]

        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[1:-1]]

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in vs))

        return [rho] + vs + [p]

    @staticmethod
    def macropri_to_macrocon(pris, cfg):
        rho, p = pris[0], pris[-1]

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[1:-1]]

        # Compute the energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in pris[1:-1])

        return [rho] + rhovs + [E]

    def set_backend(self, backend, nscalupts, nonce, linoff):
        super().set_backend(backend, nscalupts, nonce, linoff)

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.tfluxlin')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.negdivconfbgk')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.limiter')
        self._be.pointwise.register('pyfr.solvers.bgk.kernels.macrostate')

        # Setup solver matrices and parameters
        self.umat = self._be.const_matrix(self.u)
        self.Mmat = self._be.const_matrix(np.reshape(self.M, (1, -1)))
        self.niters = self.cfg.getint('solver', 'niters')

        # Get solver constants
        tau_ref = self.cfg.getfloat('constants', 'tau_ref')
        rho_ref = self.cfg.getfloat('constants', 'rho_ref')
        P_ref = self.cfg.getfloat('constants', 'P_ref')
        omega = self.cfg.getfloat('constants', 'omega')
        Pr = self.cfg.getfloat('constants', 'Pr', 1.0)
        theta_ref = P_ref/rho_ref

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims, 'nupts': self.nupts,
            'nvars': self.nvars, 'nuvars' : self.nuvars,
            'nmvars' : self.nmvars, 'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'srcex': self._src_exprs, 'pi': np.pi,
            'niters': self.niters, 'delta': self.delta,
            'tau_ref': tau_ref, 'rho_ref': rho_ref, 
            'P_ref': P_ref, 'theta_ref' : theta_ref,
            'omega' : omega, 'Pr' : Pr,
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                f=s(self.scal_upts[uin], c), F=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts'), u=self.umat
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                f=s(self._scal_qpts, c), F=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts'), u=self.umat
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                f=s(self.scal_upts[uin], l), F=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts,
                u=self.umat
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                f=s(self._scal_qpts, l), F=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts,
                u=self.umat
            )

        plocsrc = self._ploc_in_src_exprs
        plocupts = self.ploc_at('upts') if plocsrc else None
    
        self.kernels['negdivconf'] = lambda fout: self._be.kernel(
            'negdivconfbgk', tplargs=tplargs,
            dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout],
            rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, f=self._scal_upts_cpy,
            u=self.umat, M=self.Mmat
        )

        # Positivity-preserving squeeze limiter
        if self.cfg.getbool('solver', 'limiter', False) and self.basis.order != 0:
            ub = self.basis.ubasis
            tplargs['wts'] = ub.invvdm[:,0]/np.sum(ub.invvdm[:,0])

            self.kernels['limiter'] = lambda uin: self._be.kernel(
                'limiter', tplargs=tplargs,
                dims=[self.neles], f=self.scal_upts[uin]
            )
        
        # Compute and store macroscopic variables
        self.mvars = self._be.matrix((self.nupts, self.nmvars, self.neles),
                                     extent=nonce + 'mvars', tags={'align'})

        self.kernels['macrostate'] = lambda uin: self._be.kernel(
            'macrostate', tplargs=tplargs,
            dims=[self.nupts, self.neles], f=self.scal_upts[uin],
            mvars=self.mvars, u=self.umat, M=self.Mmat
        )