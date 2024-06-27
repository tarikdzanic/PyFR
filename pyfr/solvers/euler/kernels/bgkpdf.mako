<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='bgkpdf' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'
              alpha='out fpdtype_t[${str(nupts+2*nfpts)}][${str(nvars)}]'
              bounds='out fpdtype_t[${str(nvars*2)}]'
              umin='out fpdtype_t'
              du='out fpdtype_t'>

        fpdtype_t d, rcpd, v[${ndims}], p, theta;
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
        rcpd = 1.0/d;
        % for j in range(ndims):
        v[${j}] = rcpd*u[${i}][${1+j}];
        % endfor
        theta = rcpd*${c['gamma']-1}*(u[${i}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{j}]', j=ndims)});
    % elif i < nupts + nfpts:
        d = uf[${i-nupts}][0];
        rcpd = 1.0/d;
        % for j in range(ndims):
        v[${j}] = rcpd*uf[${i-nupts}][${1+j}];
        % endfor
        theta = rcpd*${c['gamma']-1}*(uf[${i-nupts}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{j}]', j=ndims)});
    % else:
        d = uf2[${i-nupts-nfpts}][0];
        rcpd = 1.0/d;
        % for j in range(ndims):
        v[${j}] = rcpd*uf2[${i-nupts-nfpts}][${1+j}];
        % endfor
        theta = rcpd*${c['gamma']-1}*(uf2[${i-nupts-nfpts}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{j}]', j=ndims)});
    % endif

        alpha[${i}][0] = d*pow(${2*pi}*theta, ${-ndims/2.0});
        alpha[${i}][1] = 1.0/(2.0*theta);
        % for j in range(ndims):
        alpha[${i}][${2+j}] = v[${j}];
        avmax = fmax(avmax, abs(v[${j}]));
        % endfor

        thetamax = fmax(thetamax, theta);
    % endfor

        // Zero bounds
        % for i in range(nvars*2):
        bounds[${i}] = 0.0;
        % endfor

        // Compute integration bounds
        umin = -avmax - ${sigma}*sqrt(thetamax);
        du = 2*abs(umin)/${(nintpts-1.0)}; // Step size
</%pyfr:kernel>
