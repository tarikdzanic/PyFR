# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.baseadvec.kernels.smats'/>

<%pyfr:kernel name='tfluxlin' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              F='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              verts='in broadcast-col fpdtype_t[${str(nverts)}][${str(ndims)}]'
              upts='in broadcast-row fpdtype_t[${str(ndims)}]'>
    // Compute the S matrices
    fpdtype_t smats[${ndims}][${ndims}], djac;
    ${pyfr.expand('calc_smats_detj', 'verts', 'upts', 'smats', 'djac')};

    // Compute and transform the fluxes
    fpdtype_t ftemp[${ndims}];
    fpdtype_t u[${ndims}];
    int fidx;

    for (int i = 0; i < ${N[0]}; i++) {
        u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

        for (int j = 0; j < ${N[1]}; j++) {
            u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

            % if ndims == 2:
            fidx = i*${N[1]} + j;
            % for m in range(ndims):
            F[${m}][fidx] = ${' + '.join(f'smats[{m}][{k}]*u[{k}]*f[fidx]' for k in range(ndims))};
            % endfor
            % else:
            for (int k = 0; k < ${N[2]}; k++) {
                u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                fidx = i*${N[1]*N[2]} + j*${N[2]} + k;
                % for m in range(ndims):
                F[${m}][fidx] = ${' + '.join(f'smats[{m}][{k}]*u[{k}]*f[fidx]' for k in range(ndims))};
                % endfor
            }
            % endif
        }
    }
</%pyfr:kernel>
