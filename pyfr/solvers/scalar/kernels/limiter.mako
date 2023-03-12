<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='cost1' params='ui, uavg, alpha'>
    alpha = (ui[0] - bounds[0])/fmax(${1E-8}, abs(uavg[0] - bounds[0]) - (ui[0] - bounds[0]));
</%pyfr:macro>
<%pyfr:macro name='cost2' params='ui, uavg, alpha'>
    alpha = (bounds[1] - ui[0])/fmax(${1E-8}, abs(bounds[1] - uavg[0]) - (bounds[1] - ui[0]));
</%pyfr:macro>

<%include file='pyfr.solvers.baseadvec.kernels.limiter'/>
<%pyfr:kernel name='limiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              x='in broadcast fpdtype_t[${str(nupts)}][${str(ndims)}]'
              bounds='inout fpdtype_t[2]'>

    ${pyfr.expand('optimize_and_limit_1', 'u' ,'x')};
    ${pyfr.expand('optimize_and_limit_2', 'u' ,'x')};
</%pyfr:kernel>
