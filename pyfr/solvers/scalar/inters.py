from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)
import math

class ScalarIntersMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.system = self.cfg.get('solver', 'system')
        if self.system == 'advection':
            subs = self.cfg.items('constants')
            subs |= dict(x='ploc[0]', y='ploc[1]', z='ploc[2]')
            subs |= dict(abs='fabs', pi=math.pi)
            vx = self.cfg.getexpr('solver', 'vx', subs=subs)
            vy = self.cfg.getexpr('solver', 'vy', subs=subs)
            if self.ndims == 3:
                vz = self.cfg.getexpr('solver', 'vz', subs=subs)
                self.v = [vx, vy, vz]
            else:
                self.v = [vx, vy]
            ploc_in_vel = any('ploc' in vv for vv in self.v)
            assert len(self.v) == self.ndims
        elif self.system == 'burgers':
            self.v = self.cfg.getliteral('solver', 'v')
            assert len(self.v) == self.ndims
        elif self.system == 'kpp':
            assert self.ndims == 2
            self.v = None
        else:
            raise ValueError(f'Unknown system: {self.system}')

class ScalarIntInters(ScalarIntersMixin, BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.intcflux')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.intcbounds')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.c, system=self.system, v=self.v)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs,
            ploc=self._ploc_at_fpts
        )

        self.kernels['comm_bounds'] = lambda: self._be.kernel(
            'intcbounds', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs
        )

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.intcboundsl')
        self.kernels['comm_bounds_l'] = lambda: self._be.kernel(
            'intcboundsl', tplargs={}, dims=[self.ninters],
            bounds_l_lhs=self._bounds_l_lhs, bounds_l_rhs=self._bounds_l_rhs
        )
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.intcboundsh')
        self.kernels['comm_bounds_h'] = lambda: self._be.kernel(
            'intcboundsh', tplargs={}, dims=[self.ninters],
            bounds_h_lhs=self._bounds_h_lhs, bounds_h_rhs=self._bounds_h_rhs
        )
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.intcboundse')
        self.kernels['comm_bounds_e'] = lambda: self._be.kernel(
            'intcboundse', tplargs={}, dims=[self.ninters],
            bounds_e_lhs=self._bounds_e_lhs, bounds_e_rhs=self._bounds_e_rhs
        )


class ScalarMPIInters(ScalarIntersMixin, BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.mpicflux')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.mpicbounds')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.c, system=self.system, v=self.v)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs,
            ploc=self._ploc_at_fpts
        )

        self.kernels['comm_bounds'] = lambda: self._be.kernel(
            'mpicbounds', tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.mpicboundsl')
        self.kernels['comm_bounds_l'] = lambda: self._be.kernel(
            'mpicboundsl', tplargs={}, dims=[self.ninters],
            bounds_l_lhs=self._bounds_l_lhs, bounds_l_rhs=self._bounds_l_rhs
        )
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.mpicboundsh')
        self.kernels['comm_bounds_h'] = lambda: self._be.kernel(
            'mpicboundsh', tplargs={}, dims=[self.ninters],
            bounds_h_lhs=self._bounds_h_lhs, bounds_h_rhs=self._bounds_h_rhs
        )
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.mpicboundse')
        self.kernels['comm_bounds_e'] = lambda: self._be.kernel(
            'mpicboundse', tplargs={}, dims=[self.ninters],
            bounds_e_lhs=self._bounds_e_lhs, bounds_e_rhs=self._bounds_e_rhs
        )


class ScalarBaseBCInters(ScalarIntersMixin, BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.bccflux')
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.bccbounds')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars, c=self.c,
                       bctype=self.type, ninters=self.ninters,
                       system=self.system, v=self.v)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs, nl=self._pnorm_lhs,
            ploc=self._ploc_at_fpts, **self._external_vals
        )

        self.kernels['comm_bounds'] = lambda: self._be.kernel(
            'bccbounds', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.bccboundsl')
        self.kernels['comm_bounds_l'] = lambda: self._be.kernel(
            'bccboundsl', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, bounds_l_lhs=self._bounds_l_lhs,
            nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
        )
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.bccboundsh')
        self.kernels['comm_bounds_h'] = lambda: self._be.kernel(
            'bccboundsh', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, bounds_h_lhs=self._bounds_h_lhs,
            nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
        )
        self._be.pointwise.register('pyfr.solvers.scalar.kernels.bccboundse')
        self.kernels['comm_bounds_e'] = lambda: self._be.kernel(
            'bccboundse', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, bounds_e_lhs=self._bounds_e_lhs,
            nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
        )

class ScalarFixedBCInters(ScalarBaseBCInters):
    type = 'fixed'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['u'], lhs)


class ScalarFreeBCInters(ScalarBaseBCInters):
    type = 'free'
    cflux_state = 'ghost'
