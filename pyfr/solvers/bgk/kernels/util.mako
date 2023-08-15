# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.matrices'/>

<%pyfr:macro name='compute_moments' params='f, u, M, w'>
    % for i in range(ndims + 2):
    w[${i}] = 0.0;
    % endfor

    fpdtype_t fm;
    for (int i = 0; i < ${nuvars}; i++) {
        fm = M[0][i]*f[i];

        w[0] += fm;
        w[1] += fm*u[i][0];
        w[2] += fm*u[i][1];
        % if ndims == 2:
        w[3] += 0.5*fm*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
        % elif ndims == 3:
        w[3] += fm*u[i][2];
        w[4] += 0.5*fm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
        % endif

        % if delta:
        w[${ndims+1}] += M[0][i]*f[i + ${nuvars}];
        % endif
    }
</%pyfr:macro>

<%pyfr:macro name='con_to_pri' params='w, q'>
    q[0] = w[0];
    q[1] = w[1]/w[0];
    q[2] = w[2]/w[0];
    % if ndims == 2:
    q[3] = ${c['gamma']-1.0}*(w[3] - 0.5*q[0]*(q[1]*q[1] + q[2]*q[2]));
    % elif ndims == 3:
    q[3] = w[3]/w[0];
    q[4] = ${c['gamma']-1.0}*(w[4] - 0.5*q[0]*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]));
    % endif
</%pyfr:macro>

<%pyfr:macro name='pri_to_con' params='q, w'>
    w[0] = q[0];
    w[1] = q[0]*q[1];
    w[2] = q[0]*q[2];
    % if ndims == 2:
    w[3] = ${1.0/(c['gamma']-1.0)}*q[3] + 0.5*q[0]*(q[1]*q[1] + q[2]*q[2]);
    % elif ndims == 3:
    w[3] = q[0]*q[3];
    w[4] = ${1.0/(c['gamma']-1.0)}*q[4] + 0.5*q[0]*(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    % endif
</%pyfr:macro>

<%pyfr:macro name='compute_alpha_Gaussian' params='q, alpha'>
    fpdtype_t theta_tmp = q[${ndims+1}]/q[0];
    alpha[0] = q[0]*pow(${2*pi}*theta_tmp, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta_tmp);
    % for i in range(ndims):
    alpha[${i+2}] = q[${i+1}];
    % endfor
</%pyfr:macro>

<%pyfr:macro name='compute_alpha_ellipsoidal' params='q, T, alpha'>
    // Compute inverse temperature tensor
    fpdtype_t Tinv[${ndims}][${ndims}] = {0}, invdet;
    % if ndims == 2:
    ${pyfr.expand('compute_2x2inverse', 'T', 'Tinv', 'invdet')};
    alpha[0] = q[0]/sqrt(${(2*pi)**ndims}/invdet);
    alpha[1] = Tinv[0][0];
    alpha[2] = Tinv[0][1];
    alpha[3] = Tinv[1][1];
    alpha[4] = q[1];
    alpha[5] = q[2];
    % elif ndims == 3:
    ${pyfr.expand('compute_3x3inverse', 'T', 'Tinv', 'invdet')};
    alpha[0] = q[0]/sqrt(${(2*pi)**ndims}/invdet);
    alpha[1] = Tinv[0][0];
    alpha[2] = Tinv[0][1];
    alpha[3] = Tinv[0][2];
    alpha[4] = Tinv[1][1];
    alpha[5] = Tinv[1][2];
    alpha[6] = Tinv[2][2];
    alpha[7] = q[1];
    alpha[8] = q[2];
    alpha[9] = q[2];
    % endif
</%pyfr:macro>

<%pyfr:macro name='compute_Maxwellian_distribution' params='alpha, u, g'>
    // Compute square of pecular velocity
    % if ndims == 2:
    fpdtype_t dv2 = (u[0]-alpha[2])*(u[0]-alpha[2]) + (u[1]-alpha[3])*(u[1]-alpha[3]);
    % elif ndims == 3:
    fpdtype_t dv2 = (u[0]-alpha[2])*(u[0]-alpha[2]) + (u[1]-alpha[3])*(u[1]-alpha[3]) + (u[2]-alpha[4])*(u[2]-alpha[4]);
    % endif

    // Compute monatomic Maxwellian
    g = alpha[0]*exp(-alpha[1]*dv2);
</%pyfr:macro>

<%pyfr:macro name='compute_ellipsoidal_distribution' params='alpha, u, g'>
    // Compute pecular velocity
    % if ndims == 2:
    fpdtype_t c0 = u[0] - alpha[4];
    fpdtype_t c1 = u[1] - alpha[5];

    fpdtype_t dv2 = alpha[1]*c0*c0 + 2*alpha[2]*c0*c1 + alpha[3]*c1*c1;
    % elif ndims == 3:
    fpdtype_t c0 = u[0] - alpha[7];
    fpdtype_t c1 = u[1] - alpha[8];
    fpdtype_t c2 = u[2] - alpha[9];

    fpdtype_t dv2 = alpha[1]*c0*c0 + 2*alpha[2]*c0*c1 + 2*alpha[3]*c0*c2 + 2*alpha[4]*c1*c1 + 2*alpha[5]*c1*c2 + alpha[6]*c2*c2;
    % endif

    // Compute ellipsoidal distribution
    g = alpha[0]*exp(-0.5*dv2);
</%pyfr:macro>

<%pyfr:macro name='compute_temperature_tensor' params='T, f, g, alpha, u, M'>
    // Compute peculiar stress tensor
    for (int i = 0; i < ${nuvars}; i++) {
        fpdtype_t Mf = M[0][i]*f[i];
        fpdtype_t Mg = M[0][i]*g[i];
        % for j,k in pyfr.ndrange(ndims, ndims):
        T[${j}][${k}] += (${1.0/Pr}*Mg + ${1.0 - 1.0/Pr}*Mf)*(u[i][${j}] - alpha[${j+2}])*(u[i][${k}] - alpha[${k+2}]);
        % endfor
    }
</%pyfr:macro>

<%pyfr:macro name='iterate_DVM_BGK' params='alpha, w, u, M'>
    fpdtype_t R[${ndims+2}];
    fpdtype_t J[${ndims+2}][${ndims+2}], Jinv[${ndims+2}][${ndims+2}];
    fpdtype_t mmnts[${ndims+2}];
    fpdtype_t gm, Mgm, invdet;

    // Pre-compute theta*delta/2.0
    fpdtype_t td2 = ${0.25*delta}/alpha[1];
    
    for (int iter = 0; iter < ${niters}; iter++) {
        // Zero cost-function and Jacobian
        % for ivar in range(ndims+2):
        R[${ivar}] = 0; 
        % for jvar in range(ndims+2):
        J[${ivar}][${jvar}] = 0; 
        % endfor
        % endfor

        // Compute discrete Maxwellian
        for (int i = 0; i < ${nuvars}; i++) {
            ${pyfr.expand('compute_Maxwellian_distribution', 'alpha', 'u[i]', 'gm')};

            // Precompute moment factors
            Mgm = M[0][i]*gm;
            mmnts[0] = Mgm;
            mmnts[1] = Mgm*u[i][0];
            mmnts[2] = Mgm*u[i][1];
            % if ndims == 2:
            mmnts[3] = 0.5*Mgm*(u[i][0]*u[i][0] + u[i][1]*u[i][1]);
            % elif ndims == 3:
            mmnts[3] = Mgm*u[i][2];
            mmnts[4] = 0.5*Mgm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2]);
            % endif

            // Add internal energy effects
            % if delta:
            mmnts[${ndims+1}] += Mgm*td2;
            % endif

            // Compute Jacobian
            % for ivar in range(ndims+2):
            R[${ivar}] += mmnts[${ivar}];

            J[${ivar}][0] += mmnts[${ivar}]/alpha[0];
            % if ndims == 2:
            J[${ivar}][1] += -mmnts[${ivar}]*( (u[i][0]-alpha[2])*(u[i][0]-alpha[2])
                                             + (u[i][1]-alpha[3])*(u[i][1]-alpha[3]) );
            % elif ndims == 3:
            J[${ivar}][1] += -mmnts[${ivar}]*( (u[i][0]-alpha[2])*(u[i][0]-alpha[2])
                                             + (u[i][1]-alpha[3])*(u[i][1]-alpha[3])
                                             + (u[i][2]-alpha[4])*(u[i][2]-alpha[4]) );
            % endif
            % for dvar in range(ndims):
            J[${ivar}][${2+dvar}] += mmnts[${ivar}]*2*alpha[1]*(u[i][${dvar}] - alpha[${2+dvar}]);
            % endfor
            % endfor
        }


        // Get defect
        % for var in range(ndims+2):
        R[${var}] -= w[${var}]; 
        % endfor

        // Compute inverse Jacobian
        % if ndims == 2:
        ${pyfr.expand('compute_4x4inverse', 'J', 'Jinv', 'invdet')};
        % elif ndims == 3:
        ${pyfr.expand('compute_5x5inverse', 'J', 'Jinv', 'invdet')};
        % endif

        // Take Newton iteration
        % for var in range(ndims+2):
        alpha[${var}] = alpha[${var}] - (${' + '.join('Jinv[{var}][{i}]*R[{i}]'.format(var=var, i=i) for i in range(ndims+2))});
        % endfor
    }
</%pyfr:macro>

<%pyfr:macro name='iterate_DVM_ESBGK' params='alpha, w, u, M'>
</%pyfr:macro>



