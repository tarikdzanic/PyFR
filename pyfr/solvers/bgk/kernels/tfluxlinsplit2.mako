# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>

<%pyfr:kernel name='tfluxlinsplit2' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(nvars)}]'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'>
    // Compute the S matrices
    fpdtype_t smats[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', 'smats', 'djac')};

    // Compute and transform the fluxes
    for (int j = 0; j < ${nuvars}; j++) {
        F[j] = ${' + '.join(f'smats[2][{k}]*u[j][{k}]*f[j]' for k in range(ndims))};

        % if delta:
        F[j + ${nuvars}] = ${' + '.join(f'smats[2][{k}]*u[j][{k}]*f[j + {nuvars}]' for k in range(ndims))};
        % endif
    }
</%pyfr:kernel>
