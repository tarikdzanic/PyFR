# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:kernel name='collision' ndim='2'
              t='scalar fpdtype_t'
              f='inout fpdtype_t[${str(nvars)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nuvars)}]'
              tau='out fpdtype_t'
              alpha='out fpdtype_t[${str(navars)}]'>

    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'f', 'u', 'M', 'w')};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Compute collision time based on viscosity model
    fpdtype_t theta;
    ${pyfr.expand('compute_tau', 'q', 'theta', 'tau')};

    // Compute equilibrium distribution function via DVM and rotational temperature
    fpdtype_t theta_rot;
    % if Pr == 1:
    ${pyfr.expand('iterate_alpha_BGK', 'f', 'w', 'q', 'u', 'M', 'alpha', 'theta_rot')};
    % else:
    ${pyfr.expand('iterate_alpha_ESBGK', 'f', 'w', 'q', 'u', 'M', 'alpha', 'theta_rot')};
    % endif

    // Set source term
    fpdtype_t g;
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u[i]', 'g')};

        // Set source
        f[i] = (g - f[i])/tau;
        % if delta:
        f[i + ${nuvars}] = (theta_rot*g - f[i + ${nuvars}])/tau;
        % endif
    }

</%pyfr:kernel>
