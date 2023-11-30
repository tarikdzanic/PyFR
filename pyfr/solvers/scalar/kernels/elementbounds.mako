<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='costl' params='ui, uavg, g, bounds'>
    g = ui[0];
</%pyfr:macro>

<%pyfr:macro name='costh' params='ui, uavg, g, bounds'>
    g = -ui[0];
</%pyfr:macro>

<%pyfr:macro name='coste' params='ui, uavg, g, bounds'>
    g = -ui[0]*ui[0];
</%pyfr:macro>

<%include file='pyfr.solvers.baseadvec.kernels.limiter'/>
<%pyfr:kernel name='elementbounds' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              x='in broadcast fpdtype_t[${str(nupts)}][${str(ndims)}]'
              bounds='out fpdtype_t[3]'
              bounds_l='out fpdtype_t[${str(nfaces)}]'
              bounds_h='out fpdtype_t[${str(nfaces)}]'
              bounds_e='out fpdtype_t[${str(nfaces)}]'>

    fpdtype_t tmp, g;
    ${pyfr.expand('optimize_bounds_l', 'u', 'tmp', 'x', 'g', 'bounds')};
    bounds[0] = g;
    ${pyfr.expand('optimize_bounds_h', 'u', 'tmp', 'x', 'g', 'bounds')};
    bounds[1] = g;
    ${pyfr.expand('optimize_bounds_e', 'u', 'tmp', 'x', 'g', 'bounds')};
    bounds[2] = g;

    % for i in range(nfaces):
    bounds_l[${i}] = bounds[0];
    bounds_h[${i}] = bounds[1];
    bounds_e[${i}] = bounds[2];
    % endfor
</%pyfr:kernel>
