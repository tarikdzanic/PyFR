<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% eps = 1E-12 %>

// Density bound
<%pyfr:macro name='g1' params='u, g, bounds'>
    g = u[0] - ${gbnds[0]};
</%pyfr:macro>

<%pyfr:macro name='g2' params='u, g, bounds'>
    fpdtype_t p = ${c['gamma'] - 1}*(u[${nvars - 1}] - 0.5/u[0]*(${pyfr.dot('u[{i}]', i=(1, ndims + 1))}));
    g = p - ${gbnds[1]};
</%pyfr:macro>

<%pyfr:macro name='h_from_g' params='gu, gavg, h'>
    h = gu >= 0 ? gu/gavg : gu/(gavg - gu);
</%pyfr:macro>

<%pyfr:macro name='cost1' params='ui, uavg, h, bounds'>
    ${pyfr.expand('g1', 'ui' ,'gu', 'bounds')};
    ${pyfr.expand('g1', 'uavg' ,'gavg', 'bounds')};
    ${pyfr.expand('h_from_g', 'gu' ,'gavg', 'h')};
</%pyfr:macro>

<%pyfr:macro name='cost2' params='ui, uavg, h, bounds'>
    ${pyfr.expand('g2', 'ui' ,'gu', 'bounds')};
    ${pyfr.expand('g2', 'uavg' ,'gavg', 'bounds')};
    if (gu >= 0) {
        h = gu/gavg;
    }
    else {
        fpdtype_t rho = ui[0];
        fpdtype_t drho = uavg[0] - ui[0];

        fpdtype_t m[${ndims}], dm[${ndims}], dm2 = 0, m2 = 0, dmm = 0;
        % for i in range(ndims):
        m[${i}] = ui[${i+1}];
        dm[${i}] = uavg[${i+1}] - ui[${i+1}];
        dm2 += dm[${i}]*dm[${i}];
        m2 += m[${i}]*m[${i}];
        dmm += dm[${i}]*m[${i}];
        % endfor

        fpdtype_t E = ui[${nvars-1}];
        fpdtype_t dE = uavg[${nvars-1}]- ui[${nvars-1}];

        fpdtype_t A = drho*dE - 0.5*dm2;
        fpdtype_t B = E*drho + rho*dE - dmm - drho*${gbnds[1]/(c['gamma'] - 1)};
        fpdtype_t C = rho*E - 0.5*m2 + rho*${gbnds[1]/(c['gamma'] - 1)};

        if (abs(A) < ${eps}) {
            h = -1;
        }
        else {
            fpdtype_t a1 = (-B + sqrt(B*B - 4*A*C))/(2*A);
            fpdtype_t a2 = (-B - sqrt(B*B - 4*A*C))/(2*A);

            h = -fmax(a1, a2);
        }
    }
</%pyfr:macro>

<%include file='pyfr.solvers.baseadvec.kernels.limiter'/>
<%pyfr:kernel name='limiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              x='in broadcast fpdtype_t[${str(nupts)}][${str(ndims)}]'
              bounds='inout fpdtype_t[3]'>
    fpdtype_t gu, gavg;
    fpdtype_t uavg[${nvars}];
    ${pyfr.expand('compute_element_average', 'u', 'uavg')};

    fpdtype_t gavg1, gavg2;
    ${pyfr.expand('g1', 'uavg', 'gavg1', 'bounds')};
    ${pyfr.expand('g2', 'uavg', 'gavg2', 'bounds')};

    if (gavg1 < ${eps} || gavg2 < ${eps}) {
        % for i,j in pyfr.ndrange(nupts, nvars):
        u[${i}][${j}] = uavg[${j}];
        % endfor
    }
    else {
        ${pyfr.expand('optimize_and_limit_1', 'u', 'uavg', 'x', 'bounds')};
        ${pyfr.expand('optimize_and_limit_2', 'u', 'uavg', 'x', 'bounds')};
    }
</%pyfr:kernel>
