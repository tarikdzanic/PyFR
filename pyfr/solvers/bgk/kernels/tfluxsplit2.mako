# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='tfluxsplit2' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>

    // Compute and transform the fluxes
    for (int j = 0; j < ${nuvars}; j++) {
        F[j] = ${' + '.join(f'smats[2][{k}]*u[j][{k}]*f[j]' for k in range(ndims))};

        % if delta:
        F[j + ${nuvars}] = ${' + '.join(f'smats[2][{k}]*u[j][{k}]*f[j + {nuvars}]' for k in range(ndims))};
        % endif
    }
</%pyfr:kernel>