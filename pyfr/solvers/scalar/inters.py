from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)
class ScalarIntersMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.system = self.cfg.get('solver', 'system')
        if self.system in ['advection', 'burgers']:
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

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.c, system=self.system, v=self.v)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class ScalarMPIInters(ScalarIntersMixin, BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.mpicflux')

        system = self.cfg.get('solver', 'system')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.c, system=self.system, v=self.v)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class ScalarBaseBCInters(ScalarIntersMixin, BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.scalar.kernels.bccflux')

        system = self.cfg.get('solver', 'system')
        tplargs = dict(ndims=self.ndims, nvars=self.nvars, c=self.c,
                       bctype=self.type, ninters=self.ninters,
                       system=self.system, v=self.v)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )

class ScalarFixedBCInters(ScalarBaseBCInters):
    type = 'fixed'

    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        self.c |= self._exp_opts(['u'], lhs)


class ScalarFreeBCInters(ScalarBaseBCInters):
    type = 'free'
    cflux_state = 'ghost'
