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
        w[3] += 0.5*fm*(u[i][0]*u[i][0] + u[i][1]*u[i][1])${f' + M[0][i]*f[i + {nuvars}]' if delta else ''};
        % elif ndims == 3:
        w[3] += fm*u[i][2];
        w[4] += 0.5*fm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])${f' + M[0][i]*f[i + {nuvars}]' if delta else ''};
        % endif
    }
</%pyfr:macro>

<%pyfr:macro name='compute_extended_moments' params='f, u, M, exm'>
    % for i in range(4*ndims-2):
    exm[${i}] = 0.0;
    % endfor

    fpdtype_t fm;
    for (int i = 0; i < ${nuvars}; i++) {
        fm = M[0][i]*f[i];

        exm[0] += fm;
        exm[1] += fm*u[i][0];
        exm[2] += fm*u[i][1];
        % if ndims == 2:
        exm[3] += fm*u[i][0]*u[i][0];
        exm[4] += fm*u[i][0]*u[i][1];
        exm[5] += fm*u[i][1]*u[i][1];
        % elif ndims == 3:
        exm[3] += fm*u[i][2];
        exm[4] += fm*u[i][0]*u[i][0];
        exm[5] += fm*u[i][0]*u[i][1];
        exm[6] += fm*u[i][0]*u[i][2];
        exm[7] += fm*u[i][1]*u[i][1];
        exm[8] += fm*u[i][1]*u[i][2];
        exm[9] += fm*u[i][2]*u[i][2];
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

<%pyfr:macro name='compute_alpha_ellipsoidal' params='q, T, Tvinv, alpha'>
    // Compute inverse temperature tensor
    fpdtype_t Tv[${ndims}][${ndims}], invdet;

    fpdtype_t theta_pr = ${1.0/Pr}*q[${ndims+1}]/q[0];
    % for i,j in pyfr.ndrange(ndims, ndims):
    Tv[${i}][${j}] = T[${i}][${j}]${' + theta_pr' if i == j else ''};
    % endfor

    % if ndims == 2:
    ${pyfr.expand('compute_2x2inverse', 'Tv', 'Tvinv', 'invdet')};
    % elif ndims == 3:
    ${pyfr.expand('compute_3x3inverse', 'Tv', 'Tvinv', 'invdet')};
    % endif
    
    alpha[0] = q[0]/sqrt(${(2*pi)**ndims}/invdet);
    alpha[1] = theta_pr;
    % for i in range(ndims):
    alpha[${2 + i}] = q[${i+1}];
    % endfor
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

<%pyfr:macro name='compute_ellipsoidal_distribution' params='alpha, u, Tvinv, g'>
    // Compute pecular velocity
    % if ndims == 2:
    fpdtype_t c0 = u[0] - alpha[2];
    fpdtype_t c1 = u[1] - alpha[3];

    fpdtype_t dv2 = Tvinv[0][0]*c0*c0 + 2*Tvinv[0][1]*c0*c1 + Tvinv[1][1]*c1*c1;
    % elif ndims == 3:
    fpdtype_t c0 = u[0] - alpha[2];
    fpdtype_t c1 = u[1] - alpha[3];
    fpdtype_t c2 = u[2] - alpha[4];

    fpdtype_t dv2 = Tvinv[0][0]*c0*c0 + 2*Tvinv[0][1]*c0*c1 + 2*Tvinv[0][2]*c0*c2 + Tvinv[1][1]*c1*c1 + 2*Tvinv[1][2]*c1*c2 + Tvinv[2][2]*c2*c2;
    % endif

    // Compute ellipsoidal distribution
    g = alpha[0]*exp(-0.5*dv2);
    //printf("%f %f %f %f %f %f %f %f %f\n", alpha[0], alpha[1], alpha[2], alpha[3], Tvinv[0][0], Tvinv[0][1], Tvinv[1][1], dv2, g);
</%pyfr:macro>

<%pyfr:macro name='compute_temperature_tensor' params='T, f, q, u, M'>
    // Compute scaled temperature tensor
    for (int i = 0; i < ${nuvars}; i++) {
        fpdtype_t Mf = ${1.0 - 1.0/Pr}*M[0][i]*f[i];
        % for j,k in pyfr.ndrange(ndims, ndims):
        T[${j}][${k}] += Mf*(u[i][${j}] - q[${j+1}])*(u[i][${k}] - q[${k+1}])/q[0];
        % endfor
    }
</%pyfr:macro>

<%pyfr:macro name='compute_temperature_Jacobian' params='T, dTvinvdtheta, theta'>
    % if ndims == 2:
    /* 
    Compute d/dtheta of T + theta*I where T is symmetric
    Use M = det(T)*T^(-1), so that T^(-1) = M/det(T)

    M[0][0] = T[1][1] + theta;
    M[0][1] = -T[0][1];
    M[1][1] = T[0][0] + theta;

    dM[0][0] = 1.0;
    dM[0][1] = 0.0;
    dM[1][1]  = 1.0;
    */

    fpdtype_t det, ddet, rcpdet2;
    det = (T[0][0] + theta)*(T[1][1] + theta) - T[0][1]*T[0][1];
    ddet = 2*theta + T[0][0] + T[1][1];
    rcpdet2 = 1.0/(det*det);

    dTvinvdtheta[0] = rcpdet2*(det - (T[1][1] + theta)*ddet);
    dTvinvdtheta[1] = rcpdet2*(T[0][1]*ddet);
    dTvinvdtheta[2] = rcpdet2*(det - (T[0][0] + theta)*ddet);
    % elif ndims == 3:
    % endif
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
            mmnts[3] = 0.5*Mgm*(u[i][0]*u[i][0] + u[i][1]*u[i][1])${' + Mgm*td2' if delta else ''};
            % elif ndims == 3:
            mmnts[3] = Mgm*u[i][2];
            mmnts[4] = 0.5*Mgm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])${' + Mgm*td2' if delta else ''};
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

<%pyfr:macro name='iterate_DVM_ESBGK' params='alpha, w, T, Tvinv, u, M'>
    fpdtype_t R[${ndims+2}];
    fpdtype_t J[${ndims+2}][${ndims+2}], Jinv[${ndims+2}][${ndims+2}];
    fpdtype_t Tv[${ndims}][${ndims}];
    fpdtype_t mmnts[${ndims+2}];
    fpdtype_t gm, Mgm, invdet;

    // Pre-compute theta*delta/2.0
    fpdtype_t td2 = ${delta*Pr/2.0}*alpha[1];

    // Compute temperature tensor
    % for i,j in pyfr.ndrange(ndims, ndims):
    Tv[${i}][${j}] = T[${i}][${j}]${' + alpha[1]' if i == j else ''};
    % endfor

    // Allocate temperature tensor Jacobian
    fpdtype_t dTvinvdtheta[${3 if ndims == 2 else 6}];
   
    for (int iter = 0; iter < ${niters}; iter++) {
        // Zero cost-function and Jacobian
        % for ivar in range(ndims+2):
        R[${ivar}] = 0; 
        % for jvar in range(ndims+2):
        J[${ivar}][${jvar}] = 0; 
        % endfor
        % endfor

        // Compute Jacobian of temperature tensor w.r.t. temperature
        ${pyfr.expand('compute_temperature_Jacobian', 'T', 'dTvinvdtheta', 'alpha[1]')};

        // Compute discrete Maxwellian
        for (int i = 0; i < ${nuvars}; i++) {
            ${pyfr.expand('compute_ellipsoidal_distribution', 'alpha', 'u[i]', 'Tvinv', 'gm')};

            // Precompute moment factors
            Mgm = M[0][i]*gm;
            mmnts[0] = Mgm;
            mmnts[1] = Mgm*u[i][0];
            mmnts[2] = Mgm*u[i][1];
            % if ndims == 2:
            mmnts[3] = 0.5*Mgm*(u[i][0]*u[i][0] + u[i][1]*u[i][1])${' + Mgm*td2' if delta else ''};
            % elif ndims == 3:
            mmnts[3] = Mgm*u[i][2];
            mmnts[4] = 0.5*Mgm*(u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])${' + Mgm*td2' if delta else ''};
            % endif

            // Compute Jacobian
            % if ndims == 2:
            fpdtype_t c0 = u[i][0] - alpha[2];
            fpdtype_t c1 = u[i][1] - alpha[3];

            % for ivar in range(ndims+2):
            R[${ivar}] += mmnts[${ivar}];

            J[${ivar}][0] += mmnts[${ivar}]/alpha[0];
            J[${ivar}][1] += mmnts[${ivar}]*(-0.5*c0*c0*dTvinvdtheta[0] - c0*c1*dTvinvdtheta[1] - 0.5*c1*c1*dTvinvdtheta[2]);
            J[${ivar}][2] += mmnts[${ivar}]*(Tvinv[0][0]*c0 + Tvinv[0][1]*c1);
            J[${ivar}][3] += mmnts[${ivar}]*(Tvinv[0][1]*c0 + Tvinv[1][1]*c1);
            % endfor
            % elif ndims == 3:
            fpdtype_t c0 = u[i][0] - alpha[2];
            fpdtype_t c1 = u[i][1] - alpha[3];
            fpdtype_t c2 = u[i][2] - alpha[4];
            % endif
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

        // Update inverse temperature tensor
        % for i,j in pyfr.ndrange(ndims, ndims):
        Tv[${i}][${j}] = T[${i}][${j}]${' + alpha[1]' if i == j else ''};
        % endfor

        % if ndims == 2:
        ${pyfr.expand('compute_2x2inverse', 'Tv', 'Tvinv', 'invdet')};
        % elif ndims == 3:
        ${pyfr.expand('compute_3x3inverse', 'Tv', 'Tvinv', 'invdet')};
        % endif
    }
</%pyfr:macro>



