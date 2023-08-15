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

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha_Gaussian', 'q', 'alpha')};
    
    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM_BGK', 'alpha', 'w', 'u', 'M')};

    // Apply ES-BGK model if necessary
    % if Pr != 1:
    fpdtype_t gm[${nuvars}];
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'gm[i]')};
    }

    // Compute initial temperature tensor
    fpdtype_t T[${ndims}][${ndims}] = {0};
    ${pyfr.expand('compute_temperature_tensor', 'T', 'f', 'gm', 'alpha', 'u', 'M')}; 

    // Get alpha vector
    fpdtype_t alpha_es[${4*ndims-2}];
    ${pyfr.expand('compute_alpha_ellipsoidal', 'q', 'T', 'alpha_es')};

    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM_ESBGK', 'alpha_es', 'w', 'u', 'M')};    
    % endif

    // Compute temperature for internal DOFs if necessary
    % if delta:
    fpdtype_t theta = q[${ndims+1}]/q[0];
    % endif

    // Set RHS state
    for (int i = 0; i < ${nuvars}; i++) {
        ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'fr[i]')};

        // Apply internal energy effects
        % if delta:
        fr[i + ${nuvars}] = fr[i]*theta*${delta/2.0};
        % endif
    }
</%pyfr:macro>
