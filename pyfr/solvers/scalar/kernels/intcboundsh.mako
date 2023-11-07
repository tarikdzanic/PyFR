<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcboundsh' ndim='1'
              bounds_h_lhs='inout view fpdtype_t'
              bounds_h_rhs='inout view fpdtype_t'>

    bounds_h_lhs = bounds_h_rhs = fmin(bounds_h_lhs, bounds_h_rhs);
</%pyfr:kernel>
