<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_pressure' params='u, p'>
    fpdtype_t rcpd = 1.0/u[0];
    fpdtype_t E = u[${nvars - 1}];

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*rcpd*(${pyfr.dot('u[{i}]', i=(1, ndims + 1))}));
</%pyfr:macro>

<%pyfr:macro name='get_minima' params='u, dmin, pmin'>
    fpdtype_t d, p;
    fpdtype_t ui[${nvars}];

    dmin = ${fpdtype_max}; pmin = ${fpdtype_max};

    for (int i = 0; i < ${nupts}; i++)
    {
        % for j in range(nvars):
        ui[${j}] = u[i][${j}];
        % endfor

        d = ui[0];
        ${pyfr.expand('compute_pressure', 'ui', 'p')};
        dmin = fmin(dmin, d); pmin = fmin(pmin, p);
    }
</%pyfr:macro>

<%pyfr:kernel name='limiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'>

    fpdtype_t dmin, pmin;
    // Check if solution is within bounds
    ${pyfr.expand('get_minima', 'u', 'dmin', 'pmin')};

    // Filter if out of bounds
    if (dmin < ${d_min} || pmin < ${p_min})
    {
        // Compute local element mean
        fpdtype_t ubar[${nvars}];
        % for i in range(nvars):
        ubar[${i}] = ${' + '.join('{jx}*u[{j}][{i}]'.format(i=i, j=j, jx=jx)
                            for j, jx in enumerate(wts) if jx != 0)};
        % endfor

        fpdtype_t dbar = ubar[0];
        fpdtype_t pbar;
        ${pyfr.expand('compute_pressure', 'ubar', 'pbar')};

        fpdtype_t alpha = 1;
        fpdtype_t tmp = (dbar - ${d_min})/fmax(dbar - dmin, ${1e-12});
        alpha = min(alpha, tmp);

        tmp = (pbar - ${p_min})/fmax(pbar - pmin, ${1e-12});
        alpha = min(alpha, tmp);

        // Apply limiting
        % for i,j in pyfr.ndrange(nupts, nvars):
        u[${i}][${j}] = alpha*u[${i}][${j}] + (1 - alpha)*ubar[${j}];
        % endfor
    }
</%pyfr:kernel>
