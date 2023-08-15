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
    // Using p = ${c['theta']}*w[0];
    w[${ndims+1}] = ${1.0/(c['gamma']-1.0)}*${c['theta']}*w[0] + (0.5/w[0])*${pyfr.dot('w[{i}]', i=(1, ndims+1))};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha', 'q', 'alpha')};
    
    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM', 'alpha', 'w', 'u', 'M')};

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
        fr[i + ${nuvars}] = fr[i]*${c['theta']*delta/2.0};
        % endif
    }
</%pyfr:macro>
