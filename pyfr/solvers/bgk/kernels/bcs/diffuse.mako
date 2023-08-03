# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr' externs='ploc, t'>
    // Get LHS conserved state
    fpdtype_t wl[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'fl', 'wl')};

    // Convert to primitives
    fpdtype_t ql[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'wl', 'ql')};

    // Compute RHS state
    fpdtype_t w[${ndims+2}] = {0};
    w[0] = wl[0];
    % for i, v in enumerate('uvw'[:ndims]):
    w[${i+1}] = w[0]*${c[v]};
    % endfor
    // Using p = ${c['theta']}*w[0];
    w[${ndims+1}] = ${1.0/(c['gamma']-1.0)}*${c['theta']}*w[0] + (0.5/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims+1))};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha', 'q', 'alpha')};
    
    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM', 'alpha', 'w')};

    // Compute mass-preserving scaling factor
    fpdtype_t u[${ndims}];
    fpdtype_t un, eta1 = 0.0, eta2 = 0.0;
    int fidx;
    for (int i = 0; i < ${N[0]}; i++) {
        u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

        for (int j = 0; j < ${N[1]}; j++) {
            u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

            % if ndims == 2:
            fidx = i*${N[1]} + j;

            // Compute equilibrium distribution at fidx-th velocity point
            ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'fr[fidx]')};

            // Balance mass flux
            un = ${pyfr.dot('u[{i}]', 'nl[{i}]', i=ndims)};
            if (un > 0.0) {
                eta1 += ${M}*fl[fidx]*abs(un);
            }
            else {
                eta2 += ${M}*fr[fidx]*abs(un);
            }
            % else:
            for (int k = 0; k < ${N[2]}; k++) {
                u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                fidx = i*${N[1]*N[2]} + j*${N[2]} + k;

                // Compute equilibrium distribution at fidx-th velocity point
                ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'fr[fidx]')};

                // Balance mass flux
                un = ${pyfr.dot('u[{i}]', 'nl[{i}]', i=ndims)};
                if (un > 0.0) {
                    eta1 += ${M}*fl[fidx]*abs(un);
                }
                else {
                    eta2 += ${M}*fr[fidx]*abs(un);
                }
            }
            % endif
        }
    }

    // Set RHS state to preserve zero mass flux
    fpdtype_t eta = eta1/eta2;
    for (int fidx = 0; fidx < ${nvars}; fidx++) {
        fr[fidx] *= eta;
    }
</%pyfr:macro>
