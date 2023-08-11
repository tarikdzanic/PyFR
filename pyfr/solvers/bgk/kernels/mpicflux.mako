# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='mpicflux' ndim='1'
              fl='inout view fpdtype_t[${str(nvars)}]'
              fr='in mpi fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>
    // Perform the Riemann solve and write out the common normal fluxes
    fpdtype_t Fn, ui[${ndims}], fli, fri;
    for (int i = 0; i < ${nuvars}; i++) {
        % for j in range(ndims):
        ui[${j}] = u[i][${j}];
        % endfor

        fli = fl[i];
        fri = fr[i];
        ${pyfr.expand('rsolve', 'fli', 'fri', 'nl', 'Fn', 'u')};

        fl[i] = magnl*Fn;

        % if delta:
        fli = fl[i + ${nuvars}];
        fri = fr[i + ${nuvars}];
        ${pyfr.expand('rsolve', 'fli', 'fri', 'nl', 'Fn', 'u')};

        fl[i + ${nuvars}] = magnl*Fn;
        % endif
    }
</%pyfr:kernel>
