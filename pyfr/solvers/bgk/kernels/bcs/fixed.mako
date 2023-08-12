# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
    fpdtype_t w[${ndims + 2}];
    w[0] = ${c['rho']};
% for i, v in enumerate('uvw'[:ndims]):
    w[${i + 1}] = (${c['rho']})*(${c[v]});
% endfor
    w[${ndims+1}] = ${c['p']}/${c['gamma'] - 1}
                    + (0.5/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims + 1))};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha', 'q', 'alpha')};

    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM', 'alpha', 'w', 'u', 'M')};

    // Compute temperature
    % if delta:
    fpdtype_t theta = q[${ndims+1}]/q[0];
    % endif

    // Precompute necessary data for Shakov model (for Prandtl number effects)
    % if Pr != 1.0:
    fpdtype_t p = q[${ndims+1}];
    fpdtype_t S[${ndims}] = {0};
    fpdtype_t Pr = ${Pr};
    ${pyfr.expand('compute_Shakov_heatflux', 'alpha', 'f', 'M', 'u', 'S')};
    % endif

    // Set RHS state
    for (int i = 0; i < ${nuvars}; i++) {
        ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'i', 'fr[i]')};

        // Apply Shakov model
        % if Pr != 1.0:
        ${pyfr.expand('apply_Shakov_model', 'alpha', 'u', 'S', 'p', 'theta', 'Pr', 'i', 'fr[i]')};
        % endif

        // Apply internal energy effects
        % if delta:
        fr[i + ${nuvars}] = fr[i]*theta*${delta/2.0};
        % endif
    }
</%pyfr:macro>
