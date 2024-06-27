<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='bgkbounds' ndim='1'
              alpha='in fpdtype_t[${str(nupts+2*nfpts)}][${str(nvars)}]'
              umin='in fpdtype_t'
              du='in fpdtype_t'
              bounds='inout fpdtype_t[${str(nvars*2)}]'>

        // Compute trapezoid rule integration weight (for interior points)
        fpdtype_t M = pow(du, ${float(ndims)});

        fpdtype_t f, g, f_min, f_max, g_min, g_max, theta; 
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
                theta = 0.5/(alpha[${i}][1]);
                g = theta*${delta/2.0}*f;

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
</%pyfr:kernel>
