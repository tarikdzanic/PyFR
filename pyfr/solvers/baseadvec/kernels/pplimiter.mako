<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% eps = 1e-12 %>
<%pyfr:macro name='compute_pressure' params='u, p'>
    fpdtype_t rcpd = 1.0/u[0];
    fpdtype_t E = u[${nvars - 1}];
    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*rcpd*(${pyfr.dot('u[{i}]', i=(1, ndims + 1))}));
</%pyfr:macro>


<%pyfr:kernel name='pplimiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>

    // Compute local element mean
    fpdtype_t ubar[${nvars}], ui[${nvars}];
    % for i in range(nvars):
    ubar[${i}] = ${' + '.join('{jx}*u[{j}][{i}]'.format(i=i, j=j, jx=jx)
                        for j, jx in enumerate(meanwts) if jx != 0)};
    % endfor

    fpdtype_t alpha = 0.0, p;
    for (int i = 0; i < ${nupts}; i++){
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor
        fpdtype_t rho = u[i][0];
        ${pyfr.expand('compute_pressure', 'ui', 'p')};

        if (p < ${p_min}) {
            fpdtype_t drho = ubar[0] - u[i][0];
            fpdtype_t m[${ndims}], dm[${ndims}], dm2 = 0, m2 = 0, dmm = 0;
            % for j in range(ndims):
            m[${j}] = u[i][${j+1}];
            dm[${j}] = ubar[${j+1}] - u[i][${j+1}];
            dm2 += dm[${j}]*dm[${j}];
            m2 += m[${j}]*m[${j}];
            dmm += dm[${j}]*m[${j}];
            % endfor
            fpdtype_t E = u[i][${nvars-1}];
            fpdtype_t dE = ubar[${nvars-1}]- u[i][${nvars-1}];
            fpdtype_t A = drho*dE - 0.5*dm2;
            fpdtype_t B = E*drho + rho*dE - dmm - drho*${p_min/(c['gamma'] - 1)};
            fpdtype_t C = rho*E - 0.5*m2 + rho*${p_min/(c['gamma'] - 1)};

            fpdtype_t tmp;
            if (abs(A) < ${eps}) {
                tmp = 1.0;
            }
            else {
                fpdtype_t a1 = (-B + sqrt(B*B - 4*A*C))/(2*A);
                fpdtype_t a2 = (-B - sqrt(B*B - 4*A*C))/(2*A);
                tmp = max(a1, a2);
            }

            alpha = max(tmp, alpha);
        }

        if (rho < ${d_min}) {
            alpha = max(alpha, (${d_min} - rho)/(ubar[0] - u[i][0]));
        }     
    }

    alpha = max(0.0, min(1.0, alpha));
    
    // Apply limiting
    if (alpha > 0) {
        % for i,j in pyfr.ndrange(nupts, nvars):
        u[${i}][${j}] = (1 - alpha)*u[${i}][${j}] + alpha*ubar[${j}];
        % endfor
    }
</%pyfr:kernel>