<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcmean' ndim='1'
              mean_lhs='inout view fpdtype_t[${str(nvars)}]'
              mean_rhs='inout view fpdtype_t[${str(nvars)}]'>

    fpdtype_t tmp;
    % for i in range(nvars):
    tmp = mean_rhs[${i}];
    mean_rhs[${i}] = mean_lhs[${i}];
    mean_lhs[${i}] = tmp;
    % endfor
</%pyfr:kernel>
