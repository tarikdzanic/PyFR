# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:kernel name='collision' ndim='2'
              t='scalar fpdtype_t'
              f='inout fpdtype_t[${str(nvars)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nuvars)}]'
              tau='out fpdtype_t'>

    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'f', 'u', 'M', 'w')};

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

    // Compute collision time based on viscosity model
    fpdtype_t p = q[${ndims+1}];
    fpdtype_t theta = p/q[0];
    % if viscosity_law == 'constant-tau':
    tau = ${tau_ref/Pr};
    % elif viscosity_law == 'constant-viscosity':
    tau = ${tau_ref*P_ref/Pr}/p;
    % elif viscosity_law == 'power-law':
    tau = ${tau_ref/Pr}*pow(theta/${theta_ref}, ${omega})/(p/${P_ref});
    % elif viscosity_law == 'sutherland':
    // mu = mu_ref*(T/T_ref)^1.5 * (T_ref + T_s)/(T + T_s)
    fpdtype_t theta_rat = theta/${theta_ref};
    tau = (${tau_ref*P_ref*(theta_ref + theta_s)/Pr}/p)*theta_rat*sqrt(theta_rat)/(theta + ${theta_s});
    % endif

    // Set source term
    fpdtype_t g;
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        % if Pr != 1:
        ${pyfr.expand('compute_ellipsoidal_distribution', 'alpha', 'u[i]', 'Tvinv', 'g')};
        % else:
        ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'g')};
        % endif

        // Set source
        f[i] = (g - f[i])/tau;
        % if delta:
        f[i + ${nuvars}] = (${delta/2.0}*theta*g - f[i + ${nuvars}])/tau;
        % endif
    }

</%pyfr:kernel>
