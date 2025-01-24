# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
    // Get LHS conserved state
    fpdtype_t wl[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'fl', 'u', 'M', 'wl')};

    // Compute RHS state
    fpdtype_t w[${ndims+2}];
    w[0] = ${c['p']}/${c['theta']};
    % for i in range(ndims):
    w[${i+1}] = wl[${i+1}];
    % endfor
    w[${ndims+1}] = ${c['p']}/${(c['gamma'] - 1)}
                    + (0.5/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims + 1))};
    
    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Compute equilibrium distribution function via DVM and rotational temperature
    fpdtype_t theta_rot, alpha[${navars}];
    ${pyfr.expand('iterate_alpha', 'f', 'w', 'q', 'u', 'M', 'alpha', 'theta_rot')};

    // Set RHS state
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u[i]', 'fr[i]')};

        // Apply internal energy effects
        % if delta:
        fr[i + ${nuvars}] = theta_rot*fr[i];
        % endif
    }
</%pyfr:macro>
