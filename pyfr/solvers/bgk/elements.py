# -*- coding: utf-8 -*-

from ctypes.wintypes import PSIZE
from pyfr.solvers.baseadvec import BaseAdvectionElements
import numpy as np
from math import gamma as gamma_func

def setup_BGK(cfg, ndims):
    Nx = cfg.getint('solver', 'Nx')
    Ny = cfg.getint('solver', 'Ny')
    if ndims == 3:
        Nz = cfg.getint('solver', 'Nz')

    vmax = cfg.getfloat('solver', 'vmax')
    u0 = cfg.getfloat('solver', 'u0')
    v0 = cfg.getfloat('solver', 'v0')
    if ndims == 3:
        w0 = cfg.getfloat('solver', 'w0')

    delta = cfg.getint('solver', 'delta', 0)
    if delta != 0:
        emax = cfg.getfloat('solver', 'emax')
        Ne = cfg.getint('solver', 'Ne')

    if ndims == 2:
        if delta:
            nvars = Nx*Ny*Ne
            u = np.zeros((nvars, ndims + 1))
            w = np.zeros((nvars))

            ux = np.linspace(-1, 1, Nx)*vmax + u0
            uy = np.linspace(-1, 1, Ny)*vmax + v0
            ue = np.linspace(0, 1, Ne)*emax

            [uxx, uyy, uee] = np.meshgrid(ux, uy, ue)
            u[:,0] = np.reshape(uxx, (-1))
            u[:,1] = np.reshape(uyy, (-1))
            u[:,2] = np.reshape(uee, (-1))

            wts_x = np.ones_like(ux)
            wts_x[0] = 0.5
            wts_x[-1] = 0.5
            wts_x = wts_x/np.sum(wts_x)*2*vmax

            wts_y = np.ones_like(uy)
            wts_y[0] = 0.5
            wts_y[-1] = 0.5
            wts_y = wts_y/np.sum(wts_y)*2*vmax

            wts_e = np.ones_like(ue)
            wts_e[0] = 0.5
            wts_e[-1] = 0.5
            wts_e = wts_e/np.sum(wts_e)*emax

            wts = np.outer(wts_e, np.outer(wts_x, wts_y))
            w = np.reshape(wts, (-1))
            PSint = w
        else:
            nvars = Nx*Ny
            u = np.zeros((nvars, ndims))
            w = np.zeros((nvars))

            ux = np.linspace(-1, 1, Nx)*vmax + u0
            uy = np.linspace(-1, 1, Ny)*vmax + v0

            [uxx, uyy] = np.meshgrid(ux, uy)
            u[:,0] = np.reshape(uxx, (-1))
            u[:,1] = np.reshape(uyy, (-1))

            wts_x = np.ones_like(ux)
            wts_x[0] = 0.5
            wts_x[-1] = 0.5
            wts_x = wts_x/np.sum(wts_x)*2*vmax

            wts_y = np.ones_like(uy)
            wts_y[0] = 0.5
            wts_y[-1] = 0.5
            wts_y = wts_y/np.sum(wts_y)*2*vmax

            wts = np.outer(wts_x, wts_y)
            w = np.reshape(wts, (-1))
            PSint = w

    elif ndims == 3:
        if delta:
            raise ValueError('Internal DOFs not implemented for 3D.')
        nvars = Nx*Ny*Nz
        u = np.zeros((nvars, ndims))
        w = np.zeros((nvars))

        ux = np.linspace(-1, 1, Nx)*vmax + u0
        uy = np.linspace(-1, 1, Ny)*vmax + v0
        uz = np.linspace(-1, 1, Nz)*vmax + w0

        [uxx, uyy, uzz] = np.meshgrid(ux, uy, uz)
        u[:,0] = np.reshape(uxx, (-1))
        u[:,1] = np.reshape(uyy, (-1))
        u[:,2] = np.reshape(uzz, (-1))

        wts_x = np.ones_like(ux)
        wts_x[0] = 0.5
        wts_x[-1] = 0.5
        wts_x = wts_x/np.sum(wts_x)*2*vmax

        wts_y = np.ones_like(uy)
        wts_y[0] = 0.5
        wts_y[-1] = 0.5
        wts_y = wts_y/np.sum(wts_y)*2*vmax

        wts_z = np.ones_like(uz)
        wts_z[0] = 0.5
        wts_z[-1] = 0.5
        wts_z = wts_z/np.sum(wts_z)*2*vmax

        wts = np.outer(wts_z, np.outer(wts_x, wts_y))
        w = np.reshape(wts, (-1))
        PSint = w
    
    psi = np.zeros((nvars, ndims+2))
    for i in range(nvars):
        psi[i, 0] = 1
        for j in range(ndims):
            psi[i, 1+j] = u[i, j]

        if delta != 0:
            psi[i, -1] = 0.5*np.linalg.norm(u[i,:-1])**2 + u[i,-1]
        else:
            psi[i, -1] = 0.5*np.linalg.norm(u[i,:])**2
    
    return [u, PSint, psi]
    
def iterate_DVM(Uloc, u, ndims, moments, PSint, gamma, niters, delta):
    if delta != 0:
        def compute_discrete_maxwellian(alpha):
            # Compute the macro/micro velocity defect
            dv2 =  (u[...,0] - alpha[2])**2
            dv2 += (u[...,1] - alpha[3])**2
            if ndims == 3:
                dv2 += (u[...,2] - alpha[4])**2

            Mv = (alpha[0]*np.exp(-alpha[1]*dv2))

            theta = 1.0/(2.0*alpha[1])
            zeta = u[...,-1]
            lam = 1.0/gamma_func(delta/2.0)
            Me = lam*(zeta/theta)**(0.5*delta - 1.)*(1./theta)*np.exp(-zeta/theta)

            M = Mv*Me

            return M
    else:
        def compute_discrete_maxwellian(alpha):
            # Compute the macro/micro velocity defect
            dv2 =  (u[...,0] - alpha[2])**2
            dv2 += (u[...,1] - alpha[3])**2
            if ndims == 3:
                dv2 += (u[...,2] - alpha[4])**2

            M = (alpha[0]*np.exp(-alpha[1]*dv2))
            return M
    
    # Change local variables into alpha vector
    if ndims == 2:
        [rholoc, rhouloc, rhovloc, Eloc] = Uloc
        p = (gamma - 1.)*(Eloc - 0.5*(rhouloc**2 + rhovloc**2)/rholoc)
        theta = p/rholoc
        alpha = np.zeros(ndims+2)

        alpha[0] = rholoc/(2*np.pi*theta)**(ndims/2.0) # A
        alpha[1] = 1.0/(2*theta) # B
        alpha[2] = rhouloc/rholoc # C,D,E
        alpha[3] = rhovloc/rholoc # C,D,E
    elif ndims == 3:
        [rholoc, rhouloc, rhovloc, rhowloc, Eloc] = Uloc
        p = (gamma - 1.)*(Eloc - 0.5*(rhouloc**2 + rhovloc**2 + rhowloc**2)/rholoc)
        theta = p/rholoc
        alpha = np.zeros(ndims+2)

        alpha[0] = rholoc/(2*np.pi*theta)**(ndims/2.0) # A
        alpha[1] = 1.0/(2*theta) # B
        alpha[2] = rhouloc/rholoc # C,D,E
        alpha[3] = rhovloc/rholoc # C,D,E
        alpha[4] = rhowloc/rholoc # C,D,E

    Mloc = compute_discrete_maxwellian(alpha)

    # Perform Newton iterations to find optimal Maxwellian
    F = [0.0]*(ndims+2)
    for _ in range(niters):
        # Derivatives with respect to alpha
        Q = [None]*(ndims+2)
        Q[0] = 1.0/alpha[0]
        if ndims == 2:
            Q[1] = -((u[...,0] - alpha[2])**2 + (u[...,1] - alpha[3])**2)
            Q[2] = 2*alpha[1]*(u[...,0] - alpha[2])
            Q[3] = 2*alpha[1]*(u[...,1] - alpha[3])
        elif ndims == 3:
            Q[1] = -((u[...,0] - alpha[2])**2 + (u[...,1] - alpha[3])**2 + (u[...,2] - alpha[4])**2)
            Q[2] = 2*alpha[1]*(u[...,0] - alpha[2])
            Q[3] = 2*alpha[1]*(u[...,1] - alpha[3])
            Q[4] = 2*alpha[1]*(u[...,2] - alpha[4])
        
        if delta:
            Q[1] += (delta - 4*u[:,-1]*alpha[1])/(2*alpha[1])
        
        J = np.zeros((ndims+2, ndims+2))
        for ivar in range(ndims+2):
            psiM = moments[:,ivar]*Mloc
            F[ivar] = np.dot(PSint, psiM) - Uloc[ivar]
            for jvar in range(ndims+2):
                J[ivar, jvar] = np.dot(PSint, Q[jvar]*psiM)

        alpha = alpha - np.linalg.inv(J) @ F
        Mloc = compute_discrete_maxwellian(alpha)
    if np.max(F) > 1e-6:
        print('Did not converge: ', F)
    return Mloc

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

        [self.u, self.PSint, self.moments] = setup_BGK(cfg, self.ndims)
        self.nvars = len(self.u)
        
        self.iterate_ICs = cfg.getbool('solver', 'iterate_ICs', True)
        self.delta = cfg.getint('solver', 'delta')

        super().__init__(basiscls, eles, cfg)

    # Initial Maxwellian state
    def pri_to_con(self, pris, cfg):
        rho, U, p = pris[0], pris[1:-1], pris[-1]
        
        # Multiply velocity components by rho
        rhoUs = [rho*c for c in U]

        # Compute the internal/total energy
        gamma = cfg.getfloat('constants', 'gamma')
        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in U)
        M = np.zeros((self.nupts, self.nvars, self.neles))
        
        # (nupts, _, nelems) = np.shape(M)

        niters = cfg.getint('solver', 'niters') if self.iterate_ICs else 0 # Large iteration count for ICs
        for uidx in range(self.nupts):
            for eidx in range(self.neles):
                # Get local variables
                rholoc = rho if np.isscalar(rho) else rho[uidx, eidx] 
                rhouloc = rhoUs[0] if np.isscalar(rhoUs[0]) else rhoUs[0][uidx, eidx] 
                rhovloc = rhoUs[1] if np.isscalar(rhoUs[1]) else rhoUs[1][uidx, eidx]
                if self.ndims == 3:
                    rhowloc = rhoUs[2] if np.isscalar(rhoUs[2]) else rhoUs[2][uidx, eidx]
                Eloc = E if np.isscalar(E) else E[uidx, eidx]

                if self.ndims == 2:
                    Uloc = [rholoc, rhouloc, rhovloc, Eloc]
                elif self.ndims == 3:
                    Uloc = [rholoc, rhouloc, rhovloc, rhowloc, Eloc]

                # Compute local Maxwellian
                M[uidx, :, eidx] = iterate_DVM(Uloc, self.u, self.ndims, self.moments, self.PSint, gamma, niters, self.delta)

        return M

    @staticmethod
    def con_to_pri(cons, cfg, PSint, u, ndims):
        f = cons

        rho = np.dot(PSint, f.swapaxes(0,1))

        # Compute velocities
        Vs = []
        for i in range(ndims):
            rhoU = np.dot(PSint, (f.swapaxes(0,2)*u[:,i]).swapaxes(1,2)).T
            Vs.append(rhoU/rho)
        
        # Compute strains
        if ndims == 2:
            sxy = np.dot(PSint, (f.swapaxes(0,2)*u[:,0]*u[:,1]).swapaxes(1,2)).T
            Ss = [sxy]
        elif ndims == 3:
            sxy = np.dot(PSint, (f.swapaxes(0,2)*u[:,0]*u[:,1]).swapaxes(1,2)).T
            sxz = np.dot(PSint, (f.swapaxes(0,2)*u[:,0]*u[:,2]).swapaxes(1,2)).T
            syz = np.dot(PSint, (f.swapaxes(0,2)*u[:,1]*u[:,2]).swapaxes(1,2)).T
            Ss = [sxy, sxz, syz]

        idofs = len(u.T) != ndims
        if idofs:
            E = np.dot(PSint, 0.5*(f.swapaxes(0,2)*np.linalg.norm(u[:,:-1], axis=1)**2).swapaxes(1,2)).T
            E += np.dot(PSint, (f.swapaxes(0,2)*u[:,-1]).swapaxes(1,2)).T
        else:
            E = np.dot(PSint, 0.5*(f.swapaxes(0,2)*np.linalg.norm(u, axis=1)**2).swapaxes(1,2)).T

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in Vs))

        return [rho] + Vs + [p]
    
    @staticmethod
    def con_to_vis(cons, cfg, PSint, u, ndims):
        f = cons

        rho = np.dot(PSint, f.swapaxes(0,1))

        # Compute velocities
        Vs = []
        for i in range(ndims):
            rhoU = np.dot(PSint, (f.swapaxes(0,2)*u[:,i]).swapaxes(1,2)).T
            Vs.append(rhoU/rho)
        
        # Compute strains
        if ndims == 2:
            sxy = np.dot(PSint, (f.swapaxes(0,2)*u[:,0]*u[:,1]).swapaxes(1,2)).T
            Ss = [sxy]
        elif ndims == 3:
            sxy = np.dot(PSint, (f.swapaxes(0,2)*u[:,0]*u[:,1]).swapaxes(1,2)).T
            sxz = np.dot(PSint, (f.swapaxes(0,2)*u[:,0]*u[:,2]).swapaxes(1,2)).T
            syz = np.dot(PSint, (f.swapaxes(0,2)*u[:,1]*u[:,2]).swapaxes(1,2)).T
            Ss = [sxy, sxz, syz]

        idofs = len(u.T) != ndims
        if idofs:
            E = np.dot(PSint, 0.5*(f.swapaxes(0,2)*np.linalg.norm(u[:,:-1], axis=1)**2).swapaxes(1,2)).T
            E += np.dot(PSint, (f.swapaxes(0,2)*u[:,-1]).swapaxes(1,2)).T
        else:
            E = np.dot(PSint, 0.5*(f.swapaxes(0,2)*np.linalg.norm(u, axis=1)**2).swapaxes(1,2)).T

        # Compute the pressure
        gamma = cfg.getfloat('constants', 'gamma')
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in Vs))

        return [rho] + Vs + [p] + Ss

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

        ub = self.basis.ubasis
        meanweights = ub.invvdm[:,0]/np.sum(ub.invvdm[:,0])
        self.niters = self.cfg.getint('solver', 'niters')

        self.umat = self._be.const_matrix(self.u)
        self.M = self._be.const_matrix(np.reshape(self.PSint, (1, -1)))
        lam = 1.0/gamma_func(self.delta/2.0) if self.delta else 1.0

        tau_ref = self.cfg.getfloat('constants', 'tau_ref')
        rho_ref = self.cfg.getfloat('constants', 'rho_ref')
        P_ref = self.cfg.getfloat('constants', 'P_ref')
        omega = self.cfg.getfloat('constants', 'omega')
        Pr = self.cfg.getfloat('constants', 'Pr', 1.0)
        theta_ref = P_ref/rho_ref

        self.nmvars = self.ndims + 2
        
        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims, 'nupts': self.nupts, 
            'nvars': self.nvars, 'nverts': len(self.basis.linspts), 
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs, 
            'u': self.u, 'moments': self.moments, 'PSint': self.PSint,
            'srcex': self._src_exprs, 'pi': np.pi,
            'niters': self.niters, 'wts': meanweights, 
            'delta': self.delta, 'lam': lam,
            'tau_ref': tau_ref, 'rho_ref': rho_ref, 
            'P_ref': P_ref, 'theta_ref' : theta_ref,
            'omega' : omega, 'Pr' : Pr,
            'nmvars' : self.nmvars
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
            u=self.umat, M=self.M
        )

        # Positivity-preserving squeeze limiter
        if self.cfg.getbool('solver', 'limiter', False):
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
            mvars=self.mvars, u=self.umat, M=self.M
        )