# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='tflux' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>
    // Compute and transform the fluxes
for (int j = 0; j < ${nuvars}; j++) {
    % for i in range(ndims):
    F[${i}][j] = ${' + '.join(f'smats[{i}][{k}]*u[j][{k}]*f[j]' for k in range(ndims))};
    % endfor

    % if delta:
    % for i in range(ndims):
    F[${i}][j + ${nuvars}] = ${' + '.join(f'smats[{i}][{k}]*u[j][{k}]*f[j + {nuvars}]' for k in range(ndims))};
    % endfor
    % endif
}

</%pyfr:kernel>
