<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% eps = 1E-8 %>
<% niters = 10 %>
<% gamma = 1 %>
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

<%pyfr:macro name='optimize' params='u, uavg, x, alpha'>
    // Create monomial basis
    fpdtype_t um[${nupts}][${nvars}];
    % for i,k in pyfr.ndrange(nupts, nvars):
    um[${i}][${k}] = ${' + '.join(f'{jx}*u[{j}][{k}]' 
                                  for j, jx in enumerate(invvdm[i]) if jx != 0)};
    % endfor

    // Find discrete minima
    fpdtype_t ui[${nvars}], xmin[${ndims}], alpha2;

    % for i in range(nvars):
    ui[${i}] = u[0][${i}];
    % endfor

    % for i in range(ndims):
    xmin[${i}] = x[0][${i}];
    % endfor

    !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha']

    for (int i = 1; i < ${nupts}; i++) {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha2']

        if (alpha2 < alpha) {
            alpha = alpha2;
            % for j in range(ndims):
            xmin[${j}] = x[i][${j}];
            % endfor
        }
    }

    // Optimize cost function to find minimum in element
    fpdtype_t ui2[${nvars}], x2[${ndims}];
    fpdtype_t J[${ndims}];
    for (int iter = 0; iter < ${niters}; iter++) {
        ${pyfr.expand('eval_monomial', 'um', 'xmin', 'ui')};
        !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha']

        // Numerically compute Jacobian
        % for i in range(ndims):
        % for j in range(ndims):
        x2[${j}] = xmin[${j}];
        % endfor

        x2[${i}] += ${eps};
        ${pyfr.expand('eval_monomial', 'um', 'x2', 'ui2')};
        !! CALL_COSTFUNCTION ['ui2', 'uavg', 'alpha2']

        J[${i}] = (alpha2 - alpha)/${eps};
        % endfor

        // Take gradient descent step
        % for i in range(ndims):
        xmin[${i}] -= ${gamma}*J[${i}];
        % endfor

        // Limit to element bounds
        % for i in range(ndims):
        xmin[${i}] = fmin(1.0, fmax(xmin[${i}], -1.0));
        % endfor
    }

    // Get location and cost function at minima
    ${pyfr.expand('eval_monomial', 'um', 'xmin', 'ui')};
    !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha']
</%pyfr:macro>

<%pyfr:macro name='optimize_and_limit' params='u, x'>
    // Compute mean modes
    fpdtype_t uavg[${nvars}];
    % for i in range(nvars):
    uavg[${i}] = ${' + '.join(f'{jx}*u[{j}][{i}]' 
                              for j, jx in enumerate(meanwts) if jx != 0)};
    % endfor

    // Optimize function
    fpdtype_t alpha;
    ${pyfr.expand('optimize', 'u', 'uavg', 'x', 'alpha')};

    // Limit
    alpha = fmin(1.0, fmax(0, -alpha));
    % for i,j in pyfr.ndrange(nupts, nvars):
    u[${i}][${j}] = (1 - alpha)*u[${i}][${j}] + alpha*uavg[${j}];
    % endfor
    
</%pyfr:macro>
