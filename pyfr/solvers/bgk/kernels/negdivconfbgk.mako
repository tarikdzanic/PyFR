# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:kernel name='negdivconfbgk' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              f='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'f', 'w')};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha', 'q', 'alpha')};

    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM', 'alpha', 'w')};

    // Compute collision time based on viscosity model
    fpdtype_t p = q[${ndims+1}];
    fpdtype_t theta = p/q[0];
    fpdtype_t tau = ${tau_ref}*pow(theta/${theta_ref}, ${omega})/(p/${P_ref});

    // Compute the BGK source term
    fpdtype_t u[${ndims}], g;
    int fidx;

    for (int i = 0; i < ${N[0]}; i++) {
        u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

        for (int j = 0; j < ${N[1]}; j++) {
            u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

            % if ndims == 2:
            fidx = i*${N[1]} + j;

            // Compute equilibrium distribution at fidx-th velocity point
            ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'g')};

            // Set source
            tdivtconf[i] = -rcpdjac*tdivtconf[i] + (g - f[i])/tau;
            % else:
            for (int k = 0; k < ${N[2]}; k++) {
                u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                fidx = i*${N[1]*N[2]} + j*${N[2]} + k;
                // Compute equilibrium distribution at fidx-th velocity point

                ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'g')};
                
                // Set source
                tdivtconf[i] = -rcpdjac*tdivtconf[i] + (g - f[i])/tau;
            }
            % endif
        }
    }
</%pyfr:kernel>
