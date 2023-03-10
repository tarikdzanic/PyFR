<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='cost' params='ui, uavg, f'>
    f = (ui[0] - ${gbnds[0]})/fmax(${1E-8}, abs(uavg[0] - ${gbnds[0]}) - (ui[0] - ${gbnds[0]}));
</%pyfr:macro>

<%include file='pyfr.solvers.baseadvec.kernels.limiter'/>
<%pyfr:kernel name='limiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              x='in broadcast fpdtype_t[${str(nupts)}][${str(ndims)}]'>

    ${pyfr.expand('optimize_and_limit', 'u' ,'x')};
</%pyfr:kernel>
