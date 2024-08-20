# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:kernel name='gtau' ndim='2'
              t='scalar fpdtype_t'
              f='in fpdtype_t[${str(nvars)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nuvars)}]'
              g='out fpdtype_t[${str(nvars)}]'
              tau='out fpdtype_t'
              alpha='out fpdtype_t[${str(ndims+2)}]'
              mvars='out fpdtype_t[${str(nmvars)}]'>

    // Navier-Stokes conserved variables
    fpdtype_t w[${ndims+2}] = {0};
    ${pyfr.expand('compute_moments', 'f', 'u', 'M', 'w')};

    // Convert to primitives
    fpdtype_t q[${ndims+2}] = {0};
    ${pyfr.expand('con_to_pri', 'w', 'q')};

    // Compute equilibrium distribution function via DVM
    % if Pr != 1:
    fpdtype_t Tvinv[${ndims}][${ndims}] = {{0}};
    ${pyfr.expand('iterate_alpha_ESBGK', 'w', 'q', 'u', 'M', 'Tvinv', 'alpha')};
    % else:
    ${pyfr.expand('iterate_alpha_BGK', 'w', 'q', 'u', 'M', 'alpha')};
    % endif

    // Compute collision time based on viscosity model
    fpdtype_t theta;
    ${pyfr.expand('compute_tau', 'q', 'theta', 'tau')};

    // Compute macroscopic state for collection
    % for i in range(ndims + 2):
    mvars[${i}] = w[${i}];
    % endfor
    % for i in range(ndims + 2, nmvars):
    mvars[${i}] = 0.0;
    % endfor

    // Set source term
    fpdtype_t gi, df, df2;
    for (int i = 0; i < ${nuvars}; i++) {
        // Compute equilibrium distribution at i-th velocity point
        % if Pr != 1:
        ${pyfr.expand('compute_ellipsoidal_distribution', 'alpha', 'u[i]', 'Tvinv', 'gi')};
        % else:
        ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'gi')};
        % endif

        // Set source
        g[i] = gi;
        % if delta:
        g[i + ${nuvars}] = (${delta/2.0}*theta*gi);
        % endif

        // H1, H2
        mvars[${ndims+2  }] += M[0][i]*f[i]*log(fmax(1E-15, f[i]));
        mvars[${ndims+2+1}] += M[0][i]*f[i + ${nuvars}]*log(fmax(1E-15, f[i + ${nuvars}]));
        df = (gi - f[i]);
        df2 = (g[i + ${nuvars}] - f[i + ${nuvars}]);
        // dH1, dH2
        mvars[${ndims+2+2}] += M[0][i]*df*log(fmax(1E-15, f[i]));
        mvars[${ndims+2+3}] += M[0][i]*df2*log(fmax(1E-15, f[i + ${nuvars}]));
        // adf1, df1**2
        mvars[${ndims+2+4}] += M[0][i]*abs(df);
        mvars[${ndims+2+5}] += M[0][i]*df*df;
    }


</%pyfr:kernel>
