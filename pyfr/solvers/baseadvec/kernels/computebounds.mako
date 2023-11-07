<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% eps = 1E-4 %>
<% niters = 3 %>
<% gamma = 1 %>
<%pyfr:macro name='eval_face_monomial' params='um, x, ui'>
    % if ndims == 2:
    % for i in range(nvars):
    ui[${i}] = 0.0;
    % for j in range(nfptsperface):
    ui[${i}] += um[${j}][${i}]*pow(x, ${float(j)});
    % endfor 
    % endfor 
    % endif
</%pyfr:macro>

<%pyfr:macro name='optimize_face' params='u, x, alpha'>
    % if ndims == 2:
    // Create monomial face basis
    fpdtype_t um[${nfptsperface}][${nvars}];
    % for i,k in pyfr.ndrange(nfptsperface, nvars):
    um[${i}][${k}] = ${' + '.join(f'{jx}*u[{j}][{k}]' 
                                  for j, jx in enumerate(faceinvvdm[i]) if jx != 0)};
    % endfor

    // Find discrete minima
    fpdtype_t ui[${nvars}], xmin, alpha2;
    % for i in range(nvars):
    ui[${i}] = u[0][${i}];
    % endfor
    xmin = x[0][0];

    !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha']

    for (int i = 1; i < ${nfptsperface}; i++) {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha2']

        if (alpha2 < alpha) {
            alpha = alpha2;
            xmin = x[i][0];
        }
    }

    // Optimize cost function to find minimum in element
    fpdtype_t ui2[${nvars}], x2, dx;
    fpdtype_t J, H;
    fpdtype_t dalpha[3];

    for (int iter = 0; iter < ${niters}; iter++) {
        ${pyfr.expand('eval_face_monomial', 'um', 'xmin', 'ui')};
        !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha']

        // Numerically compute Jacobian/Hessian
        % for i in range(3):
        x2 = xmin + ${eps*(i-1)};
        ${pyfr.expand('eval_face_monomial', 'um', 'x2', 'ui2')};
        !! CALL_COSTFUNCTION ['ui2', 'uavg', 'alpha2']
        dalpha[${i}] = alpha2;
        % endfor
        J = (dalpha[2] - dalpha[0])/${2*eps};
        H = (dalpha[2] - 2*dalpha[1] + dalpha[0])/${eps**2};

        // Compute Newton search direction
        dx = J/H;

        // Take Newton step
        xmin -= dx;

        // Limit to element bounds
        xmin = fmin(1.0, fmax(xmin, -1.0));
    }

    // Get location and cost function at minima
    ${pyfr.expand('eval_face_monomial', 'um', 'xmin', 'ui')};
    !! CALL_COSTFUNCTION ['ui', 'uavg', 'alpha']
    % endif
</%pyfr:macro>
