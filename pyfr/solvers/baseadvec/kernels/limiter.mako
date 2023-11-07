<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% dxi = 1E-4 %>
<% eps = 1E-12 %>
<% beta = 0.1 %>

<%pyfr:macro name='compute_element_average' params='u, uavg'>
    // Compute mean modes
    % for i in range(nvars):
    uavg[${i}] = ${' + '.join(f'{jx}*u[{j}][{i}]' 
                              for j, jx in enumerate(meanwts) if jx != 0)};
    % endfor
</%pyfr:macro>

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
    % if element_type == 'quad':
    // Project step from x0 -> x1 along vector dx = x1 - x0 onto a quad (-1, 1)^2
    fpdtype_t a1 = 1.0, a2 = 1.0;

    if (dx[0] > ${eps} && x1[0] > 1.0) {
        a1 = (1.0 - x0[0])/dx[0];
    }
    else if (dx[0] < ${-eps} && x1[0] < -1.0) {
        a1 = (-1.0 - x0[0])/dx[0];
    }

    if (dx[1] > ${eps} && x1[1] > 1.0) {
        a2 = (1.0 - x0[1])/dx[1];
    }
    else if (dx[1] < ${-eps} && x1[1] < -1.0) {
        a2 = (-1.0 - x0[1])/dx[1];
    }

    fpdtype_t adx = max(0.0, min(a1, a2));
    % for i in range(ndims):
    x1[${i}] = fmax(-1.0, fmin(1.0, x0[${i}] + adx*dx[${i}]));
    % endfor
    % elif element_type == 'tri':
    // Project step from x0 -> x1 along vector dx = x1 - x0 onto a triangle defined by (-1, -1), (1, -1), (-1, 1)
    fpdtype_t a1 = 1.0, a2 = 1.0, a3 = 1.0;
    if (dx[0] < ${-eps} && x1[0] < -1.0) {
        a1 = (-1.0 - x0[0])/dx[0];
    }
    if (dx[1] < ${-eps} && x1[1] < -1.0) {
        a2 = (-1.0 - x0[1])/dx[1];
    }
    // Diagonal is defined by y = -x, step will intersect if (x0, y0) in triangle and dy > -dx
    if (dx[1] > -dx[0] - ${eps} && x1[1] > -x1[0] - ${eps}) {
        // Intersection defined by solving x = -y -> x0 + a*dx = -(y0 + a*dy) -> a = -(y0 + x0)/(dx + dy)
        a3 = -(x0[0] + x0[1])/fmax(${eps}, dx[0] + dx[1]);
    }

    // Compute interestion and explicitly enforce x > -1, y > -1, y < -x
    fpdtype_t adx = max(0.0, min(a3, min(a1, a2)));
    % for i in range(ndims):
    x1[${i}] = fmax(-1.0, fmin(1.0, x0[${i}] + adx*dx[${i}]));
    % endfor
    x1[1] = min(x1[1], -x1[0]);
    % endif
</%pyfr:macro>

<%pyfr:macro name='optimize' params='u, uavg, x, hstar, bounds'>
    // Create monomial basis
    fpdtype_t um[${nupts}][${nvars}];
    % for i,k in pyfr.ndrange(nupts, nvars):
    um[${i}][${k}] = ${' + '.join(f'{jx}*u[{j}][{k}]' 
                                  for j, jx in enumerate(invvdm[i]) if jx != 0)};
    % endfor

    // Find discrete minima
    fpdtype_t ui[${nvars}], xmin[${ndims}], xmin_old[${ndims}], hmax, hmin, h2;

    % for i in range(nvars):
    ui[${i}] = u[0][${i}];
    % endfor

    % for i in range(ndims):
    xmin[${i}] = x[0][${i}];
    % endfor

    !! CALL_COSTFUNCTION ['ui', 'uavg', 'hstar', 'bounds']
    hmax = hstar;
    hmin = hstar;

    for (int i = 1; i < ${nupts}; i++) {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        !! CALL_COSTFUNCTION ['ui', 'uavg', 'h2', 'bounds']

        if (h2 < hstar) {
            hstar = h2;
            % for j in range(ndims):
            xmin[${j}] = x[i][${j}];
            % endfor
        }
        if (h2 > hmax) {
            hmax = h2;
        }
    }

    if ((hmax - hstar) < ${eps}){
        hstar = -1;
    }
    % if niters > 0:
    else {
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
            !! CALL_COSTFUNCTION ['ui', 'uavg', 'hstar', 'bounds']
        
            // Numerically compute Jacobian/Hessian
            // Compute perturbations 
            % for i,j in pyfr.ndrange(3, 3):
            x2[0] = xmin[0] + ${dxi*(i-1)};
            x2[1] = xmin[1] + ${dxi*(j-1)};
            ${pyfr.expand('eval_monomial', 'um', 'x2', 'ui2')};
            !! CALL_COSTFUNCTION ['ui2', 'uavg', 'h2', 'bounds']
            dh[${i}][${j}] = h2;
            % endfor

            // Compute Jacobian
            J[0] = (dh[2][1] - dh[0][1])/${2*dxi};
            J[1] = (dh[1][2] - dh[1][0])/${2*dxi};

            // Compute gradient descent step
            % for i in range(ndims):
            dx1[${i}] = ${-beta}*J[${i}];
            % endfor

            // Zero outward facing component if on element boundary
            % if element_type == 'quad':
            // Zero positive dx if on right boundary
            dx1[0] = xmin[0] > ${1-eps}  ? min(dx1[0], 0.0) : dx1[0];
            // Zero negative dx if on left boundary
            dx1[0] = xmin[0] < ${-1+eps} ? max(dx1[0], 0.0) : dx1[0];
            // Zero positive dy if on top boundary
            dx1[1] = xmin[1] > ${1-eps}  ? min(dx1[1], 0.0) : dx1[1];
            // Zero negative dy if on bottom boundary
            dx1[1] = xmin[1] < ${-1+eps} ? max(dx1[1], 0.0) : dx1[1];
            % elif element_type == 'tri':
            // Zero negative dx if on left boundary
            dx1[0] = xmin[0] < ${-1+eps} ? max(dx1[0], 0.0) : dx1[0];
            // Zero negative dy if on bottom boundary
            dx1[1] = xmin[1] < ${-1+eps} ? max(dx1[1], 0.0) : dx1[1];
            // Zero outward normal (1,1) if on diagonal
            if (xmin[1] > -xmin[0] - ${eps} && dx1[0] + dx1[1] > 0) {
                dx1[0] = ${1 - 2**0.5}*dx1[0];
                dx1[1] = ${1 - 2**0.5}*dx1[1];
            }
            % endif

            // Take step and project to element bounds
            % for i in range(ndims):
            xmin1[${i}] = xmin[${i}] - dx1[${i}];
            % endfor
            ${pyfr.expand('project_step_to_element', 'xmin', 'xmin1', 'dx1')};
    
            // Evaluate solution at new point
            ${pyfr.expand('eval_monomial', 'um', 'xmin1', 'ui2')};
            !! CALL_COSTFUNCTION ['ui2', 'uavg', 'hstar', 'bounds']
            % for i in range(ndims):
            xmin[${i}] = xmin1[${i}];
            % endfor

            // Attempt Newton step and compare to GD
            H[0][0] = (dh[2][1] - 2*dh[1][1] + dh[0][1]           )/${dxi**2};
            H[1][1] = (dh[1][2] - 2*dh[1][1] + dh[1][0]           )/${dxi**2};
            H[0][1] = (dh[2][2] -   dh[0][2] - dh[2][0] + dh[0][0])/${4*dxi**2};
            H[1][0] = H[0][1];

            // Invert Hessian
            det = H[0][0]*H[1][1] - H[0][1]*H[1][0];
            if (abs(det) > ${eps}) {
                invdet = 1.0/det;
                invH[0][0] =  invdet*H[1][1];
                invH[0][1] = -invdet*H[0][1];
                invH[1][0] = -invdet*H[1][0];
                invH[1][1] =  invdet*H[0][0];

                // Take Newton and gradient descent step
                % for i in range(ndims):
                dx2[${i}] = -(${' + '.join(f'invH[{i}][{j}]*J[{j}]' for j in range(ndims))});
                % endfor

                // Zero outward facing component if on element boundary
                % if element_type == 'quad':
                // Zero positive dx if on right boundary
                dx2[0] = xmin[0] > ${1-eps}  ? min(dx2[0], 0.0) : dx2[0];
                // Zero negative dx if on left boundary
                dx2[0] = xmin[0] < ${-1+eps} ? max(dx2[0], 0.0) : dx2[0];
                // Zero positive dy if on top boundary
                dx2[1] = xmin[1] > ${1-eps}  ? min(dx2[1], 0.0) : dx2[1];
                // Zero negative dy if on bottom boundary
                dx2[1] = xmin[1] < ${-1+eps} ? max(dx2[1], 0.0) : dx2[1];
                % elif element_type == 'tri':
                // Zero negative dx if on left boundary
                dx2[0] = xmin[0] < ${-1+eps} ? max(dx2[0], 0.0) : dx2[0];
                // Zero negative dy if on bottom boundary
                dx2[1] = xmin[1] < ${-1+eps} ? max(dx2[1], 0.0) : dx2[1];
                // Zero outward normal (1,1) if on diagonal
                if (xmin[1] > -xmin[0] - ${eps} && dx2[0] + dx2[1] > 0) {
                    dx2[0] = ${1 - 2**0.5}*dx2[0];
                    dx2[1] = ${1 - 2**0.5}*dx2[1];
                }
                % endif

                // Take step and project to element bounds
                % for i in range(ndims):
                xmin2[${i}] = xmin[${i}] - dx2[${i}];
                % endfor
                ${pyfr.expand('project_step_to_element', 'xmin', 'xmin2', 'dx2')};

                ${pyfr.expand('eval_monomial', 'um', 'xmin2', 'ui2')};
                !! CALL_COSTFUNCTION ['ui2', 'uavg', 'h2', 'bounds']

                if (h2 < hstar) {
                    hstar = h2;
                    % for i in range(ndims):
                    xmin[${i}] = xmin2[${i}];
                    % endfor
                }
            }

            // Track minimum value across iterations
            hmin = fmin(hmin, hstar);
        }

        // Extrapolate lower bound for hstar
        fpdtype_t dhJ2 = (${' + '.join(f'pow(J[{i}]*J[{i}]*(xmin[{i}] - xmin_old[{i}]), 2.0)' for i in range(ndims))});
        hstar = fmax(-1, hmin - sqrt(dhJ2));
    }
    % endif
</%pyfr:macro>



<%pyfr:macro name='optimize_bounds' params='u, uavg, x, hstar, bounds'>
    // Create monomial basis
    fpdtype_t um[${nupts}][${nvars}];
    % for i,k in pyfr.ndrange(nupts, nvars):
    um[${i}][${k}] = ${' + '.join(f'{jx}*u[{j}][{k}]' 
                                  for j, jx in enumerate(invvdm[i]) if jx != 0)};
    % endfor
    // Find discrete minima
    fpdtype_t ui[${nvars}], xmin[${ndims}], xmin_old[${ndims}], hmax, hmin, h2;
    % for i in range(nvars):
    ui[${i}] = u[0][${i}];
    % endfor
    % for i in range(ndims):
    xmin[${i}] = x[0][${i}];
    % endfor
    !! CALL_COSTFUNCTION ['ui', 'uavg', 'hstar', 'bounds']
    hmin = hstar;
    for (int i = 1; i < ${nupts}; i++) {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor
        !! CALL_COSTFUNCTION ['ui', 'uavg', 'h2', 'bounds']
        if (h2 < hstar) {
            hstar = h2;
            % for j in range(ndims):
            xmin[${j}] = x[i][${j}];
            % endfor
        }
    }
    % if niters > 0:
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
            !! CALL_COSTFUNCTION ['ui', 'uavg', 'hstar', 'bounds']
        
            // Numerically compute Jacobian/Hessian
            // Compute perturbations 
            % for i,j in pyfr.ndrange(3, 3):
            x2[0] = xmin[0] + ${dxi*(i-1)};
            x2[1] = xmin[1] + ${dxi*(j-1)};
            ${pyfr.expand('eval_monomial', 'um', 'x2', 'ui2')};
            !! CALL_COSTFUNCTION ['ui2', 'uavg', 'h2', 'bounds']
            dh[${i}][${j}] = h2;
            % endfor
            // Compute Jacobian
            J[0] = (dh[2][1] - dh[0][1])/${2*dxi};
            J[1] = (dh[1][2] - dh[1][0])/${2*dxi};
            // Compute gradient descent step
            % for i in range(ndims):
            dx1[${i}] = ${-beta}*J[${i}];
            % endfor
            // Zero outward facing component if on element boundary
            % if element_type == 'quad':
            // Zero positive dx if on right boundary
            dx1[0] = xmin[0] > ${1-eps}  ? min(dx1[0], 0.0) : dx1[0];
            // Zero negative dx if on left boundary
            dx1[0] = xmin[0] < ${-1+eps} ? max(dx1[0], 0.0) : dx1[0];
            // Zero positive dy if on top boundary
            dx1[1] = xmin[1] > ${1-eps}  ? min(dx1[1], 0.0) : dx1[1];
            // Zero negative dy if on bottom boundary
            dx1[1] = xmin[1] < ${-1+eps} ? max(dx1[1], 0.0) : dx1[1];
            % elif element_type == 'tri':
            // Zero negative dx if on left boundary
            dx1[0] = xmin[0] < ${-1+eps} ? max(dx1[0], 0.0) : dx1[0];
            // Zero negative dy if on bottom boundary
            dx1[1] = xmin[1] < ${-1+eps} ? max(dx1[1], 0.0) : dx1[1];
            // Zero outward normal (1,1) if on diagonal
            if (xmin[1] > -xmin[0] - ${eps} && dx1[0] + dx1[1] > 0) {
                dx1[0] = ${1 - 2**0.5}*dx1[0];
                dx1[1] = ${1 - 2**0.5}*dx1[1];
            }
            % endif
            // Take step and project to element bounds
            % for i in range(ndims):
            xmin1[${i}] = xmin[${i}] - dx1[${i}];
            % endfor
            ${pyfr.expand('project_step_to_element', 'xmin', 'xmin1', 'dx1')};
    
            // Evaluate solution at new point
            ${pyfr.expand('eval_monomial', 'um', 'xmin1', 'ui2')};
            !! CALL_COSTFUNCTION ['ui2', 'uavg', 'hstar', 'bounds']
            % for i in range(ndims):
            xmin[${i}] = xmin1[${i}];
            % endfor
            // Attempt Newton step and compare to GD
            H[0][0] = (dh[2][1] - 2*dh[1][1] + dh[0][1]           )/${dxi**2};
            H[1][1] = (dh[1][2] - 2*dh[1][1] + dh[1][0]           )/${dxi**2};
            H[0][1] = (dh[2][2] -   dh[0][2] - dh[2][0] + dh[0][0])/${4*dxi**2};
            H[1][0] = H[0][1];
            // Invert Hessian
            det = H[0][0]*H[1][1] - H[0][1]*H[1][0];
            if (abs(det) > ${eps}) {
                invdet = 1.0/det;
                invH[0][0] =  invdet*H[1][1];
                invH[0][1] = -invdet*H[0][1];
                invH[1][0] = -invdet*H[1][0];
                invH[1][1] =  invdet*H[0][0];
                // Take Newton and gradient descent step
                % for i in range(ndims):
                dx2[${i}] = -(${' + '.join(f'invH[{i}][{j}]*J[{j}]' for j in range(ndims))});
                % endfor
                // Zero outward facing component if on element boundary
                % if element_type == 'quad':
                // Zero positive dx if on right boundary
                dx2[0] = xmin[0] > ${1-eps}  ? min(dx2[0], 0.0) : dx2[0];
                // Zero negative dx if on left boundary
                dx2[0] = xmin[0] < ${-1+eps} ? max(dx2[0], 0.0) : dx2[0];
                // Zero positive dy if on top boundary
                dx2[1] = xmin[1] > ${1-eps}  ? min(dx2[1], 0.0) : dx2[1];
                // Zero negative dy if on bottom boundary
                dx2[1] = xmin[1] < ${-1+eps} ? max(dx2[1], 0.0) : dx2[1];
                % elif element_type == 'tri':
                // Zero negative dx if on left boundary
                dx2[0] = xmin[0] < ${-1+eps} ? max(dx2[0], 0.0) : dx2[0];
                // Zero negative dy if on bottom boundary
                dx2[1] = xmin[1] < ${-1+eps} ? max(dx2[1], 0.0) : dx2[1];
                // Zero outward normal (1,1) if on diagonal
                if (xmin[1] > -xmin[0] - ${eps} && dx2[0] + dx2[1] > 0) {
                    dx2[0] = ${1 - 2**0.5}*dx2[0];
                    dx2[1] = ${1 - 2**0.5}*dx2[1];
                }
                % endif
                // Take step and project to element bounds
                % for i in range(ndims):
                xmin2[${i}] = xmin[${i}] - dx2[${i}];
                % endfor
                ${pyfr.expand('project_step_to_element', 'xmin', 'xmin2', 'dx2')};
                ${pyfr.expand('eval_monomial', 'um', 'xmin2', 'ui2')};
                !! CALL_COSTFUNCTION ['ui2', 'uavg', 'h2', 'bounds']
                if (h2 < hstar) {
                    hstar = h2;
                    % for i in range(ndims):
                    xmin[${i}] = xmin2[${i}];
                    % endfor
                }
            }
            // Track minimum value across iterations
            hmin = fmin(hmin, hstar);
        }
        hstar = hmin;
        % endif
</%pyfr:macro>

<%pyfr:macro name='optimize_and_limit' params='u, uavg, x, bounds'>
    // Optimize function
    fpdtype_t hstar;
    ${pyfr.expand('optimize', 'u', 'uavg', 'x', 'hstar', 'bounds')};

    // Limit
    fpdtype_t alpha = fmin(1.0, fmax(0, -hstar));
    % for i,j in pyfr.ndrange(nupts, nvars):
    u[${i}][${j}] = (1 - alpha)*u[${i}][${j}] + alpha*uavg[${j}];
    % endfor
</%pyfr:macro>
