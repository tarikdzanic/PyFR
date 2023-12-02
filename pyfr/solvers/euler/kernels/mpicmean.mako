<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='mpicent' ndim='1'
              mean_lhs='inout view fpdtype_t[${str(nvars)}]'
              mean_rhs='in mpi fpdtype_t[${str(nvars)}]'>
    
    % for i in range(nvars):
    mean_lhs[${i}] = mean_rhs[${i}];
    % endfor
</%pyfr:kernel>
