<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='bgkbounds' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              bounds='out fpdtype_t[${str(nvars*2)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'>

        fpdtype_t alpha[${nupts + 2*nfpts}][${nvars}];
        fpdtype_t d, v[${ndims}], p, theta[${nupts + 2*nfpts}];
        fpdtype_t avmax = 0, thetamax = 0;

        // Compute interior flux point values
        fpdtype_t uf2[${nfpts}][${nvars}];
        for (int fidx = 0; fidx < ${nfpts}; fidx++) {
            for (int vidx = 0; vidx < ${nvars}; vidx++) {
                uf2[fidx][vidx] = ${pyfr.dot('m0[fidx][{k}]', 'u[{k}][vidx]', k=nupts)};
            }
        }

        // Compute alpha vector and theta over upts+fpts
    % for i in range(nupts + 2*nfpts):
    % if i < nupts:
        d = u[${i}][0];
        % for j in range(ndims):
        v[${j}] = u[${i}][${1+j}]/d;
        % endfor
        p = ${c['gamma']-1}*(u[${i}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{i}]', i=ndims)});
    % elif i < nupts + nfpts:
        d = uf[${i}][0];
        % for j in range(ndims):
        v[${j}] = uf[${i}][${1+j}]/d;
        % endfor
        p = ${c['gamma']-1}*(uf[${i}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{i}]', i=ndims)});
    % else:
        d = uf2[${i}][0];
        % for j in range(ndims):
        v[${j}] = uf2[${i}][${1+j}]/d;
        % endfor
        p = ${c['gamma']-1}*(uf2[${i}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{i}]', i=ndims)});
    % endif
        theta[${i}] = p/d;

        alpha[${i}][0] = d*pow(${2*pi}*theta[${i}], ${-ndims/2.0});
        alpha[${i}][1] = 1.0/(2.0*theta[${i}]);
        % for j in range(ndims):
        alpha[${i}][${2+j}] = v[${j}];
        avmax = fmax(avmax, abs(v[${j}]));
        % endfor

        thetamax = fmax(thetamax, theta[${i}]);
    % endfor

        // Zero bounds
        % for i in range(nvars*2):
        bounds[${i}] = 0.0;
        % endfor

        // Compute integration bounds
        fpdtype_t umin = -avmax - ${sigma}*sqrt(thetamax);
        fpdtype_t dumm = 2*(avmax + ${sigma}*sqrt(thetamax));
        fpdtype_t du = dumm/${(nintpts-1.0)}; // Step size
      
        // Compute trapezoid rule integration weight (for interior points)
        fpdtype_t M = pow(du, ${float(ndims)});

        fpdtype_t f, g, f_min, f_max, g_min, g_max; 
        fpdtype_t uu[${ndims}], mi;
        % if ndims == 2:
        for (int i = 0; i < ${nintpts}; i++) {
            uu[0] = umin + i*du;
            for (int j = 0; j < ${nintpts}; j++) {
                uu[1] = umin + j*du;

                // Compute local integration weight
                mi = M;
                mi = (i == 0 || i == ${nintpts-1}) ? 0.5*mi : mi;
                mi = (j == 0 || j == ${nintpts-1}) ? 0.5*mi : mi;

                f_min = ${fpdtype_max}; f_max = ${-fpdtype_max};
                g_min = ${fpdtype_max}; g_max = ${-fpdtype_max};
                % for i in range(nupts + 2*nfpts):
                % if ndims == 2:
                f = alpha[${i}][0]*exp(-alpha[${i}][1]*(
                        (uu[0] - alpha[${i}][2])*(uu[0] - alpha[${i}][2]) +
                        (uu[1] - alpha[${i}][3])*(uu[1] - alpha[${i}][3])));
                % elif ndims == 3:
                f = alpha[${i}][0]*exp(-alpha[${i}][1]*(
                        (uu[0] - alpha[${i}][2])*(uu[0] - alpha[${i}][2]) +
                        (uu[1] - alpha[${i}][3])*(uu[1] - alpha[${i}][3]) +
                        (uu[2] - alpha[${i}][4])*(uu[2] - alpha[${i}][4])));
                % endif
                g = theta[${i}]*${delta/2.0}*f;

                f_min = fmin(f_min, (1.0 - ${r_fac})*f);
                f_max = fmax(f_max, (1.0 + ${r_fac})*f);

                g_min = fmin(g_min, (1.0 - ${r_fac})*g);
                g_max = fmax(g_max, (1.0 + ${r_fac})*g);
                % endfor

                // Integrate bounds (rho_min, ru_min, rv_min, E_min, rho_max, ru_max, rv_max, E_max)
                bounds[0] += f_min*mi; bounds[4] += f_max*mi;
                for (int d = 0; d < ${ndims}; d++) {
                    if (uu[d] < 0) {
                        bounds[1+d] += f_max*mi*uu[d]; bounds[5+d] += f_min*mi*uu[d];
                    }
                    else {
                        bounds[1+d] += f_min*mi*uu[d]; bounds[5+d] += f_max*mi*uu[d];
                    }
                }
                bounds[3] += f_min*mi*0.5*(uu[0]*uu[0] + uu[1]*uu[1]); bounds[7] += f_max*mi*0.5*(uu[0]*uu[0] + uu[1]*uu[1]);
                bounds[3] += g_min*mi;                                 bounds[7] += g_max*mi;
            }
        }
        % endif

        // Apply minimum density bound
        bounds[0] = fmax(bounds[0], ${d_min});
</%pyfr:kernel>
