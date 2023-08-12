# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:kernel name='negdivconfbgk' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              f='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nuvars)}]'>

    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'f', 'u', 'M', 'w')};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Get alpha vector
    fpdtype_t alpha[${ndims+2}];
    ${pyfr.expand('compute_alpha', 'q', 'alpha')};

    // Compute discretely conservative equilibrium state
    ${pyfr.expand('iterate_DVM', 'alpha', 'w', 'u', 'M')};

    // Compute collision time based on viscosity model
    fpdtype_t p = q[${ndims+1}];
    fpdtype_t theta = p/q[0];
    fpdtype_t tau = ${tau_ref}*pow(theta/${theta_ref}, ${omega})/(p/${P_ref});

    // Set source term
    fpdtype_t g;
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'i', 'g')};

        // Set source
        tdivtconf[i] = -rcpdjac*tdivtconf[i] + (g - f[i])/tau;
        % if delta:
        tdivtconf[i + ${nuvars}] = -rcpdjac*tdivtconf[i + ${nuvars}] + (${delta/2.0}*theta*g - f[i + ${nuvars}])/tau;
        % endif
    }

</%pyfr:kernel>
