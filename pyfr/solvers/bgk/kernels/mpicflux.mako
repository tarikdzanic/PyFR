# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='mpicflux' ndim='1'
              fl='inout view fpdtype_t[${str(nvars)}]'
              fr='in mpi fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'>
    // Compute local velocity
    fpdtype_t u[${ndims}], Fn, fli, fri;
    int fidx;

    for (int i = 0; i < ${N[0]}; i++) {
        u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

        for (int j = 0; j < ${N[1]}; j++) {
            u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

            % if ndims == 2:
            fidx = i*${N[1]} + j;

            // Compute the Riemann solve
            fli = fl[fidx];
            fri = fr[fidx];

            ${pyfr.expand('rsolve', 'fli', 'fri', 'nl', 'Fn', 'u')};

            // Scale and write out the common normal fluxes
            fl[fidx] =  magnl*Fn;
            % else:
            for (int k = 0; k < ${N[2]}; k++) {
                u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                // Compute the Riemann solve
                fli = fl[fidx];
                fri = fr[fidx];

                ${pyfr.expand('rsolve', 'fli', 'fri', 'nl', 'Fn', 'u')};

                // Scale and write out the common normal fluxes
                fl[fidx] =  magnl*Fn;
            }
            % endif
        }
    }
</%pyfr:kernel>
