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
    for (int i = threadIdx.y; i < ${nuvars}; i += blockDim.y) {
        fpdtype_t un = ${pyfr.dot('u[i][{j}]', 'nl[{j}]', j=ndims)};
        fpdtype_t mun = magnl*un;

        if (un > 0.0) {
            fl[i] = mun*fl[i];
            % if delta:
            fl[i + ${nuvars}] = mun*fl[i + ${nuvars}];
            % endif
        }
        else {
            fl[i] = mun*fr[i];
            % if delta:
            fl[i + ${nuvars}] = mun*fr[i + ${nuvars}];
            % endif
        }
    }
</%pyfr:kernel>
