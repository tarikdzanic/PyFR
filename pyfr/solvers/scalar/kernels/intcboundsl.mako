<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcboundsl' ndim='1'
              bounds_l_lhs='inout view fpdtype_t'
              bounds_l_rhs='inout view fpdtype_t'>

    bounds_l_lhs = bounds_l_rhs = fmin(bounds_l_lhs, bounds_l_rhs);
</%pyfr:kernel>
