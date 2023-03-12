<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='costl' params='ui, uavg, alpha'>
    alpha = ui[0];
</%pyfr:macro>

<%pyfr:macro name='costh' params='ui, uavg, alpha'>
    alpha = -ui[0];
</%pyfr:macro>

<%include file='pyfr.solvers.baseadvec.kernels.limiter'/>
<%pyfr:kernel name='elementbounds' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              x='in broadcast fpdtype_t[${str(nupts)}][${str(ndims)}]'
              bounds='out fpdtype_t[2]'
              bounds_l='out fpdtype_t[${str(nfaces)}]'
              bounds_h='out fpdtype_t[${str(nfaces)}]'>

    fpdtype_t tmp, ulow, uhigh;
    ${pyfr.expand('optimize_l', 'u', 'tmp', 'x', 'ulow')};
    bounds[0] = ulow;
    ${pyfr.expand('optimize_h', 'u', 'tmp', 'x', 'uhigh')};
    bounds[1] = -uhigh;

    % if not face_bounds:
    % for i in range(nfaces):
    bounds_l[${i}] = bounds[0];
    bounds_h[${i}] = bounds[1];
    % endfor
    % endif
</%pyfr:kernel>
