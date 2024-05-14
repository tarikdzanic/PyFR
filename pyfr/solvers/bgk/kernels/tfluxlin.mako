# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>

<%pyfr:kernel name='tfluxlin' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>
    // Compute the S matrices
    fpdtype_t smats[${ndims}][${ndims}], djac, smu[${ndims}];
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', 'smats', 'djac')};

    // Compute and transform the fluxes
    for (int j = 0; j < ${nuvars}; j++) {
        % for i in range(ndims):
        smu[${i}] = ${' + '.join(f'smats[{i}][{k}]*u[j][{k}]' for k in range(ndims))};
        F[${i}][j] = smu[${i}]*f[j];
        % endfor

        % if delta:
        % for i in range(ndims):
        F[${i}][j + ${nuvars}] = smu[${i}]*f[j + ${nuvars}];
        % endfor
        % endif
    }
</%pyfr:kernel>
