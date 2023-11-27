<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% eps = 1E-12 %>

<%pyfr:macro name='g1' params='u, g, bounds'>
    % if global_bounds:
    g = u[0] - ${gbnds[0]};
    % else:
    g = u[0] - bounds[0];
    % endif
</%pyfr:macro>

<%pyfr:macro name='g2' params='u, g, bounds'>
    % if global_bounds:
    g = ${gbnds[1]} - u[0];
    % else:
    g = -bounds[1] - u[0];
    % endif
</%pyfr:macro>

<%pyfr:macro name='g3' params='u, g, bounds'>
    % if global_bounds:
    g = 0;
    % else:
    g = bounds[2] - (-0.5*u[0]*u[0]); // bounds[2] is negative since finding min(-u*u)
    % endif
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
    ${pyfr.expand('h_from_g', 'gu' ,'gavg', 'h')};
</%pyfr:macro>

<%pyfr:macro name='cost3' params='ui, uavg, h, bounds'>
    ${pyfr.expand('g3', 'ui' ,'gu', 'bounds')};
    ${pyfr.expand('g3', 'uavg' ,'gavg', 'bounds')};
    ${pyfr.expand('h_from_g', 'gu' ,'gavg', 'h')};
</%pyfr:macro>

<%include file='pyfr.solvers.baseadvec.kernels.limiter'/>
<%pyfr:kernel name='limiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              x='in broadcast fpdtype_t[${str(nupts)}][${str(ndims)}]'
              bounds='inout fpdtype_t[3]'>
    fpdtype_t gu, gavg;
    fpdtype_t uavg[${nvars}];
    ${pyfr.expand('compute_element_average', 'u', 'uavg')};

    % if apply_entropy_bounds:
    fpdtype_t gavg1, gavg2, gavg3;
    ${pyfr.expand('g1', 'uavg', 'gavg1', 'bounds')};
    ${pyfr.expand('g2', 'uavg', 'gavg2', 'bounds')};
    ${pyfr.expand('g3', 'uavg', 'gavg3', 'bounds')};

    if (gavg1 < ${eps} || gavg2 < ${eps} || gavg3 < ${eps}) {
        % for i,j in pyfr.ndrange(nupts, nvars):
        u[${i}][${j}] = uavg[${j}];
        % endfor
    }
    else {
        ${pyfr.expand('optimize_and_limit_1', 'u', 'uavg', 'x', 'bounds')};
        ${pyfr.expand('optimize_and_limit_2', 'u', 'uavg', 'x', 'bounds')};
        ${pyfr.expand('optimize_and_limit_3', 'u', 'uavg', 'x', 'bounds')};
    }
    % else:
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
    % endif
</%pyfr:kernel>
