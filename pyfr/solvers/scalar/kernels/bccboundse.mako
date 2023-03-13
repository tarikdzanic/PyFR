<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.scalar.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccboundse' ndim='1'
              ul='in view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              bounds_e_lhs='inout view fpdtype_t'>
</%pyfr:kernel>
