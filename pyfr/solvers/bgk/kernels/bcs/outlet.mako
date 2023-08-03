# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr' externs='ploc, t'>
    // Get LHS conserved state
    fpdtype_t wl[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'fl', 'wl')};

    // Compute RHS state
    fpdtype_t w[${ndims+2}];
    w[0] = wl[0];
% for i in range(ndims):
    w[${i+1}] = wl[${i+1}];
% endfor
    w[${ndims+1}] = ${c['p']}/${(c['gamma'] - 1)}
                    + (0.5/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims + 1))};
    
    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha', 'q', 'alpha')};
    
    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM', 'alpha', 'w')};

    // Set RHS state
    fpdtype_t u[${ndims}];
    int fidx;
    for (int i = 0; i < ${N[0]}; i++) {
        u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

        for (int j = 0; j < ${N[1]}; j++) {
            u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

            % if ndims == 2:
            fidx = i*${N[1]} + j;

            // Compute equilibrium distribution at fidx-th velocity point
            ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'fr[fidx]')};

            % else:
            for (int k = 0; k < ${N[2]}; k++) {
                u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                fidx = i*${N[1]*N[2]} + j*${N[2]} + k;

                // Compute equilibrium distribution at fidx-th velocity point
                ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'fr[fidx]')};
            }
            % endif
        }
    }
</%pyfr:macro>
