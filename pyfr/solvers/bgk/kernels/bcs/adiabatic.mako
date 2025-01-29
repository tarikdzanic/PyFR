# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>
    // Get LHS conserved state
    fpdtype_t wl[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'fl', 'u', 'M', 'wl')};

    // Convert to primitives
    fpdtype_t ql[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'wl', 'ql')};

    // Compute RHS state
    fpdtype_t w[${ndims+2}] = {0};
    w[0] = wl[0];
    % for i, v in enumerate('uvw'[:ndims]):
    w[${i+1}] = w[0]*${c[v]};
    % endfor
    w[${ndims+1}] = ${1.0/(c['gamma']-1.0)}*ql[${ndims+1}] + (0.5/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims+1))};

    // Get wall temperature if using interal energy model
    fpdtype_t theta = ql[${ndims+1}]/ql[0];

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Compute Maxwellian distribution function via DVM and rotational temperature
    fpdtype_t theta_rot, alpha[${ndims+2}];
    ${pyfr.expand('iterate_alpha_BGK', 'f', 'w', 'q', 'u', 'M', 'alpha', 'theta_rot')};

    // Compute mass-preserving scaling factor
    fpdtype_t Mw[${nuvars}];
    fpdtype_t un, eta1 = 0.0, eta2 = 0.0;
    for (int i = 0; i < ${nuvars}; i++) {
        un = ${pyfr.dot('u[i][{j}]', 'nl[{j}]', j=ndims)};

        // Compute equilibrium distribution at i-th velocity point
        ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'Mw[i]')};

        // Balance mass flux
        if (un > 0.0) {
            eta1 += fl[i]*M[0][i]*abs(un);
        }
        else {
            eta2 += Mw[i]*M[0][i]*abs(un);
        }
    }

    // Scale RHS state to preserve zero mass flux
    fpdtype_t eta = eta1/eta2;
    for (int i = 0; i < ${nuvars}; i++) {
        fr[i] = eta*Mw[i];
        
        // Apply internal energy effects
        % if delta:
        fr[i + ${nuvars}] = theta_rot*fr[i];
        % endif
    }
</%pyfr:macro>
