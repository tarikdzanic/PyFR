<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='mpicbounds' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='in mpi fpdtype_t[${str(nvars)}]'>

    // Swap interface solutions
    % for i in range(nvars):
    ul[${i}] = ur[${i}];
    % endfor
</%pyfr:kernel>
