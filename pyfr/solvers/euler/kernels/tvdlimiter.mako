<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='tvdlimiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ubar_int='in fpdtype_t[${str(nfaces)}][${str(nvars)}]'>


    // Compute local element mean
    fpdtype_t ubar[${nvars}];
    % for i in range(nvars):
    ubar[${i}] = ${' + '.join('{jx}*u[{j}][{i}]'.format(i=i, j=j, jx=jx)
                        for j, jx in enumerate(wts) if jx != 0)};
    % endfor

    // Compute lower/upper bounds based on min/max of means across interfaces
    fpdtype_t ul[${nvars}], uh[${nvars}];
    % for i in range(nvars):
    ul[${i}] = uh[${i}] = ubar[${i}];
    % for j in range(nfaces):
    ul[${i}] = fmin(ul[${i}], ubar_int[${j}][${i}]);
    uh[${i}] = fmax(uh[${i}], ubar_int[${j}][${i}]);
    % endfor
    % endfor

    // Find minima/maxima in element
    fpdtype_t umin[${nvars}], umax[${nvars}];
    % for i in range(nvars):
    umin[${i}] =  ${10**10};
    umax[${i}] = ${-10**10};
    % endfor

    % for i,j in pyfr.ndrange(nupts, nvars):
    umin[${j}] = fmin(umin[${j}], u[${i}][${j}]);
    umax[${j}] = fmax(umax[${j}], u[${i}][${j}]);
    % endfor


    ## // Find limiting factor
    ## fpdtype_t beta = 1;
    ## fpdtype_t tmp;
    ## % for i in range(nvars):
    ## tmp = (uh[${i}] - ubar[${i}])/fmax(umax[${i}] - ubar[${i}], ${1e-12});
    ## beta = fmin(beta, tmp);
    ## tmp = (ubar[${i}] - ul[${i}])/fmax(ubar[${i}] - umin[${i}], ${1e-12});
    ## beta = fmin(beta, tmp);
    ## % endfor
    ## beta = fmax(0, beta);

    ## // Apply limiting
    ## % for i,j in pyfr.ndrange(nupts, nvars):
    ## u[${i}][${j}] = beta*u[${i}][${j}] + (1 - beta)*ubar[${j}];
    ## % endfor

    
    fpdtype_t beta;
    fpdtype_t tmp;
    % for i in range(nvars):
    beta = 1;
    tmp = fmax(${0}, uh[${i}] - ubar[${i}])/fmax(umax[${i}] - ubar[${i}], ${1e-12});
    beta = fmin(beta, tmp);
    tmp = fmax(${0}, ubar[${i}] - ul[${i}])/fmax(ubar[${i}] - umin[${i}], ${1e-12});
    beta = fmin(beta, tmp);
    beta = fmax(0, beta);

    // Apply limiting
    % for j in range(nupts):
    u[${j}][${i}] = beta*u[${j}][${i}] + (1 - beta)*ubar[${i}];
    % endfor
    % endfor
    
    ## printf("%f %f %f %f %f %f %f %f %f %f\n", ubar[0], ul[0], uh[0], beta, 
    ## (uh[${0}] - ubar[${0}])/fmax(umax[${0}] - ubar[${0}], ${1e-12}), 
    ## (ubar[${0}] - ul[${0}])/fmax(ubar[${0}] - umin[${0}], ${1e-12}),
    ## uh[0], umax[0], ul[0], umin[0]);
</%pyfr:kernel>
