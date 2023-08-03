# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.matrices'/>

<%pyfr:macro name='compute_moments' params='f, w'>
    % for i in range(ndims + 2):
    w[${i}] = 0.0;
    % endfor

    fpdtype_t u[${ndims}], fm;
    int fidx;
    for (int i = 0; i < ${N[0]}; i++) {
        u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

        for (int j = 0; j < ${N[1]}; j++) {
            u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

            % if ndims == 2:
            fidx = i*${N[1]} + j;
            fm = ${M}*f[fidx];
            w[0] += fm;
            w[1] += fm*u[0];
            w[2] += fm*u[1];
            w[3] += 0.5*fm*(u[0]*u[0] + u[1]*u[1]);
            % else:
            for (int k = 0; k < ${N[2]}; k++) {
                u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                fidx = i*${N[1]*N[2]} + j*${N[2]} + k;
                fm = ${M}*f[fidx];
                w[0] += fm;
                w[1] += fm*u[0];
                w[2] += fm*u[1];
                w[3] += fm*u[2];
                w[4] += 0.5*fm*(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
            }
            % endif
        }
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

<%pyfr:macro name='compute_alpha' params='q, alpha'>
    fpdtype_t theta_tmp = q[${ndims+1}]/q[0];
    alpha[0] = q[0]*pow(${2*pi}*theta_tmp, ${-ndims/2.0});
    alpha[1] = 1.0/(2.0*theta_tmp);
    % for i in range(ndims):
    alpha[${i+2}] = q[${i+1}];
    % endfor
</%pyfr:macro>

<%pyfr:macro name='compute_equilibrium_distribution' params='alpha, ui, g'>
    // Compute square of pecular velocity
    % if ndims == 2:
    fpdtype_t dv2 = (ui[0]-alpha[2])*(ui[0]-alpha[2]) + (ui[1]-alpha[3])*(ui[1]-alpha[3]);
    % elif ndims == 3:
    fpdtype_t dv2 = (ui[0]-alpha[2])*(ui[0]-alpha[2]) + (ui[1]-alpha[3])*(ui[1]-alpha[3]) + (ui[2]-alpha[4])*(ui[2]-alpha[4]);
    % endif

    // Compute monatomic Maxwellian
    g = alpha[0]*exp(-alpha[1]*dv2);
</%pyfr:macro>

<%pyfr:macro name='iterate_DVM' params='alpha, w'>
    fpdtype_t R[${ndims+2}];
    fpdtype_t J[${ndims+2}][${ndims+2}], Jinv[${ndims+2}][${ndims+2}];
    fpdtype_t mmnts[${ndims+2}], u[${ndims}];
    fpdtype_t gm;
    int fidx;
    
    for (int iter = 0; iter < ${niters}; iter++) {
        // Zero cost-function and Jacobian
    % for ivar in range(ndims+2):
        R[${ivar}] = 0; 
        % for jvar in range(ndims+2):
        J[${ivar}][${jvar}] = 0; 
        % endfor
    % endfor

        for (int i = 0; i < ${N[0]}; i++) {
            u[0] = ${ubounds[0][0]} + ${(ubounds[0][1] - ubounds[0][0])/(N[0] - 1)}*i;

            for (int j = 0; j < ${N[1]}; j++) {
                u[1] = ${ubounds[1][0]} + ${(ubounds[1][1] - ubounds[1][0])/(N[1] - 1)}*j;

                % if ndims == 2:
                fidx = i*${N[1]} + j;
                ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'gm')};

                // Precompute moment factors
                fpdtype_t Mgm = ${M}*gm;
                mmnts[0] = Mgm;
                mmnts[1] = Mgm*u[0];
                mmnts[2] = Mgm*u[1];
                mmnts[3] = 0.5*Mgm*(u[0]*u[0] + u[1]*u[1]);
  
                % for ivar in range(ndims+2):
                R[${ivar}] += mmnts[${ivar}];

                J[${ivar}][0] += mmnts[${ivar}]/alpha[0];
                J[${ivar}][1] += -mmnts[${ivar}]*( (u[0]-alpha[2])*(u[0]-alpha[2])
                                                 + (u[1]-alpha[3])*(u[1]-alpha[3]) );
                J[${ivar}][2] += mmnts[${ivar}]*2*alpha[1]*(u[0] - alpha[2]);
                J[${ivar}][3] += mmnts[${ivar}]*2*alpha[1]*(u[1] - alpha[3]);
                % endfor
                % else:
                for (int k = 0; k < ${N[2]}; k++) {
                    u[2] = ${ubounds[2][0]} + ${(ubounds[2][1] - ubounds[2][0])/(N[2] - 1)}*k;

                    fidx = i*${N[1]*N[2]} + j*${N[2]} + k;
                    ${pyfr.expand('compute_equilibrium_distribution', 'alpha', 'u', 'gm')};

                    // Precompute moment factors
                    fpdtype_t Mgm = ${M}*gm;
                    mmnts[0] = Mgm;
                    mmnts[1] = Mgm*u[0];
                    mmnts[2] = Mgm*u[1];
                    mmnts[3] = Mgm*u[2];
                    mmnts[4] = 0.5*Mgm*(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
      
                    % for ivar in range(ndims+2):
                    R[${ivar}] += mmnts[${ivar}];

                    J[${ivar}][0] += mmnts[${ivar}]/alpha[0];
                    J[${ivar}][1] += -mmnts[${ivar}]*( (u[0]-alpha[2])*(u[0]-alpha[2])
                                                     + (u[1]-alpha[3])*(u[1]-alpha[3])
                                                     + (u[2]-alpha[4])*(u[2]-alpha[4]) );
                    J[${ivar}][2] += mmnts[${ivar}]*2*alpha[1]*(u[0] - alpha[2]);
                    J[${ivar}][3] += mmnts[${ivar}]*2*alpha[1]*(u[1] - alpha[3]);
                    J[${ivar}][4] += mmnts[${ivar}]*2*alpha[1]*(u[2] - alpha[4]);
                    % endfor
                }
            % endif
            }
        }

        // Get defect
        % for var in range(ndims+2):
        R[${var}] -= w[${var}]; 
        % endfor

        // Compute inverse Jacobian
        % if ndims == 2:
        ${pyfr.expand('compute_4x4inverse', 'J', 'Jinv')};
        % elif ndims == 3:
        ${pyfr.expand('compute_5x5inverse', 'J', 'Jinv')};
        % endif

        // Take Newton iteration
        % for var in range(ndims+2):
        alpha[${var}] = alpha[${var}] - (${' + '.join('Jinv[{var}][{i}]*R[{i}]'.format(var=var, i=i) for i in range(ndims+2))});
        % endfor
    }
</%pyfr:macro>



