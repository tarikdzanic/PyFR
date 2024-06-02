# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.bgk.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccflux' ndim='1'
              fl='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              magnl='in fpdtype_t'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nuvars)}]'>
    

    // Temp commit for free BCs

    // Perform the Riemann solve and write out the common normal fluxes
    fpdtype_t Fn, ui[${ndims}], fli, fri;
    for (int i = 0; i < ${nuvars}; i++) {
        % for j in range(ndims):
        ui[${j}] = u[i][${j}];
        % endfor

        fli = fl[i];
        ${pyfr.expand('rsolve', 'fli', 'fli', 'nl', 'Fn', 'ui')};
        fl[i] = magnl*Fn;

        % if delta:
        fli = fl[i + ${nuvars}];
        ${pyfr.expand('rsolve', 'fli', 'fli', 'nl', 'Fn', 'ui')};

        fl[i + ${nuvars}] = magnl*Fn;
        % endif
    }
</%pyfr:kernel>
