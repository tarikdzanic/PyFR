# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>


<%pyfr:macro name='compute_tflux' params='f, F, smats, u, fluxdim'>
    // Compute and transform the fluxes
    for (int j = 0; j < ${nuvars}; j++) {
        F[j] = ${' + '.join(f'smats[fluxdim][{k}]*u[j][{k}]*f[j]' for k in range(ndims))};

        % if delta:
        F[j + ${nuvars}] = ${' + '.join(f'smats[fluxdim][{k}]*u[j][{k}]*f[j + {nuvars}]' for k in range(ndims))};
        % endif
    }
</%pyfr:macro>

<%pyfr:kernel name='tfluxsplit0' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>

    int fluxdim = 0;
    ${pyfr.expand('compute_tflux', 'f', 'F', 'smats', 'u', 'fluxdim')};
</%pyfr:kernel>

<%pyfr:kernel name='tfluxsplit1' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>

    int fluxdim = 1;
    ${pyfr.expand('compute_tflux', 'f', 'F', 'smats', 'u', 'fluxdim')};
</%pyfr:kernel>

<%pyfr:kernel name='tfluxsplit2' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>

    int fluxdim = 2;
    ${pyfr.expand('compute_tflux', 'f', 'F', 'smats', 'u', 'fluxdim')};
</%pyfr:kernel>
