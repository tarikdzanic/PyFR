<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% dxi = 1E-4 %>
<% eps = 1E-12 %>
<% beta = 0.1 %>
<%pyfr:macro name='eval_monomial' params='um, x, ui'>
    fpdtype_t tmp;
    % for i in range(nvars):
    ui[${i}] = 0.0;
    % for j in range(nupts):
    tmp = 1.0;
    % for k in range(ndims):
    tmp *= pow(x[${k}], ${float(mdegs[j][k])});
    % endfor
    ui[${i}] += um[${j}][${i}]*tmp;
    % endfor
    % endfor
</%pyfr:macro>

<%pyfr:macro name='project_step_to_element' params='x0, x1, dx'>
    % for i in range(ndims):
    x1[${i}] = fmax(-1.0, fmin(1.0, x1[${i}]));
    % endfor
</%pyfr:macro>

<%pyfr:macro name='optimize' params='u, uavg, x, hstar'>
    // Create monomial basis
    fpdtype_t um[${nupts}][${nvars}];
    % for i,k in pyfr.ndrange(nupts, nvars):
    um[${i}][${k}] = ${' + '.join(f'{jx}*u[{j}][{k}]' 
                                  for j, jx in enumerate(invvdm[i]) if jx != 0)};
    % endfor

    // Find discrete minima
    fpdtype_t ui[${nvars}], xmin[${ndims}], xmin_old[${ndims}], h2;

    % for i in range(nvars):
    ui[${i}] = u[0][${i}];
    % endfor

    % for i in range(ndims):
    xmin[${i}] = x[0][${i}];
    % endfor

    !! CALL_COSTFUNCTION ['ui', 'uavg', 'hstar']

    for (int i = 1; i < ${nupts}; i++) {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        !! CALL_COSTFUNCTION ['ui', 'uavg', 'h2']

        if (h2 < hstar) {
            hstar = h2;
            % for j in range(ndims):
            xmin[${j}] = x[i][${j}];
            % endfor
        }
    }

    // Optimize cost function to find minimum in element
    fpdtype_t ui2[${nvars}], x2[${ndims}];
    fpdtype_t J[${ndims}], H[${ndims}][${ndims}];
    fpdtype_t invH[${ndims}][${ndims}], det, invdet;
    fpdtype_t dx1[${ndims}], dx2[${ndims}], xmin1[${ndims}], xmin2[${ndims}];
    % if ndims == 2:
    fpdtype_t dh[3][3];
    % elif ndims == 3:
    fpdtype_t dh[3][3][3];
    % endif

    for (int iter = 0; iter < ${niters}; iter++) {
        // Set old xmin
        % for i in range(ndims):
        xmin_old[${i}] = xmin[${i}];
        % endfor

        ${pyfr.expand('eval_monomial', 'um', 'xmin', 'ui')};
        !! CALL_COSTFUNCTION ['ui', 'uavg', 'hstar']
    
        // Numerically compute Jacobian/Hessian
        // Compute perturbations 
        % for i,j in pyfr.ndrange(3, 3):
        x2[0] = xmin[0] + ${dxi*(i-1)};
        x2[1] = xmin[1] + ${dxi*(j-1)};
        ${pyfr.expand('eval_monomial', 'um', 'x2', 'ui2')};
        !! CALL_COSTFUNCTION ['ui2', 'uavg', 'h2']
        dh[${i}][${j}] = h2;
        % endfor

        // Compute Jacobian
        J[0] = (dh[2][1] - dh[0][1])/${2*dxi};
        J[1] = (dh[1][2] - dh[1][0])/${2*dxi};

            ////// ADD BRANCH HERE FOR IF ON BOUNDARY OR NOT

        // Take gradient descent step
        % for i in range(ndims):
        dx1[${i}] = ${-beta}*J[${i}];
        % endfor
        // Take step and project to element bounds
        % for i in range(ndims):
        xmin1[${i}] = xmin[${i}] - dx1[${i}];
        % endfor
        ${pyfr.expand('project_step_to_element', 'xmin', 'xmin1', 'dx1')};
        
        ${pyfr.expand('eval_monomial', 'um', 'xmin1', 'ui2')};
        !! CALL_COSTFUNCTION ['ui2', 'uavg', 'hstar']
        % for i in range(ndims):
        xmin[${i}] = xmin1[${i}];
        % endfor

        // Attempt Newton step and compare to GD
        H[0][0] = (dh[2][1] - 2*dh[1][1] + dh[0][1]           )/${dxi**2};
        H[1][1] = (dh[1][2] - 2*dh[1][1] + dh[1][0]           )/${dxi**2};
        H[0][1] = (dh[2][2] -   dh[0][2] - dh[2][0] + dh[0][0])/${4*dxi**2};
        H[1][0] = H[0][1];

        // ADD CHECK HERE FOR HESSIAN DETERMINANT
        // Invert Hessian
        det = H[0][0]*H[1][1] - H[0][1]*H[1][0];
        if (abs(det) > ${eps}) {
            invdet = 1.0/det;
            invH[0][0] =  invdet*H[1][1];
            invH[0][1] = -invdet*H[0][1];
            invH[1][0] = -invdet*H[1][0];
            invH[1][1] =  invdet*H[0][0];

                ////// ADD BRANCH HERE FOR IF ON BOUNDARY OR NOT

            // Take Newton and gradient descent step
            % for i in range(ndims):
            dx2[${i}] = -(${' + '.join(f'invH[{i}][{j}]*J[{j}]' for j in range(ndims))});
            % endfor

            // Take step and project to element bounds
            % for i in range(ndims):
            xmin2[${i}] = xmin[${i}] - dx2[${i}];
            % endfor
            ${pyfr.expand('project_step_to_element', 'xmin', 'xmin2', 'dx2')};

            ${pyfr.expand('eval_monomial', 'um', 'xmin2', 'ui2')};
            !! CALL_COSTFUNCTION ['ui2', 'uavg', 'h2']

            if (h2 < hstar) {
                hstar = h2;
                % for i in range(ndims):
                xmin[${i}] = xmin2[${i}];
                % endfor
            }
        }
    }

    // Extrapolate lower bound for hstar
    fpdtype_t dhJ = (${' + '.join(f'abs(J[{i}]*(xmin[{i}] - xmin_old[{i}]))' for i in range(ndims))});
    hstar = fmax(-1, hstar - dhJ);
</%pyfr:macro>

<%pyfr:macro name='optimize_and_limit' params='u, x'>
    // Compute mean modes
    fpdtype_t uavg[${nvars}];
    % for i in range(nvars):
    uavg[${i}] = ${' + '.join(f'{jx}*u[{j}][{i}]' 
                              for j, jx in enumerate(meanwts) if jx != 0)};
    % endfor

    // Optimize function
    fpdtype_t hstar;
    ${pyfr.expand('optimize', 'u', 'uavg', 'x', 'hstar')};

    // Limit
    fpdtype_t alpha = fmin(1.0, fmax(0, -hstar));
    % for i,j in pyfr.ndrange(nupts, nvars):
    u[${i}][${j}] = (1 - alpha)*u[${i}][${j}] + alpha*uavg[${j}];
    % endfor
</%pyfr:macro>
