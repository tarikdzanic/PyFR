<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcboundse' ndim='1'
              bounds_e_lhs='inout view fpdtype_t'
              bounds_e_rhs='inout view fpdtype_t'>

    bounds_e_lhs = bounds_e_rhs = fmin(bounds_e_lhs, bounds_e_rhs);
</%pyfr:kernel>
