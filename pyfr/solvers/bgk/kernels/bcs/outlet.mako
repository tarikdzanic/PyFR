# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
    // Get LHS conserved state
    fpdtype_t wl[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'fl', 'u', 'M', 'wl')};

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

    // Apply ES-BGK model if necessary
    % if Pr != 1:
    // Compute initial temperature tensor (scaled by 1 - 1/Pr)
    fpdtype_t T[${ndims}][${ndims}] = {{0}};
    ${pyfr.expand('compute_temperature_tensor', 'T', 'f', 'q', 'u', 'M')}; 

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}], Tvinv[${ndims}][${ndims}] = {{0}};
    ${pyfr.expand('compute_alpha_ellipsoidal', 'q', 'T', 'Tvinv', 'alpha')};

    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM_ESBGK', 'alpha', 'w', 'T', 'Tvinv', 'u', 'M')};
    % else:
    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha_Gaussian', 'q', 'alpha')};

    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM_BGK', 'alpha', 'w', 'u', 'M')};
    % endif

    // Compute temperature for internal DOFs if necessary
    % if delta:
    fpdtype_t theta = q[${ndims+1}]/q[0];
    % endif

    // Set RHS state
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        % if Pr != 1:
        ${pyfr.expand('compute_ellipsoidal_distribution', 'alpha_es', 'u[i]', 'fr[i]')};
        % else:
        ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'fr[i]')};
        % endif

        // Apply internal energy effects
        % if delta:
        fr[i + ${nuvars}] = fr[i]*theta*${delta/2.0};
        % endif
    }
</%pyfr:macro>
