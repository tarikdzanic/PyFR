<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='intcbounds' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'>

    // Swap interface solutions
    fpdtype_t tmp;
    % for i in range(nvars):
    tmp = ul[${i}];
    ul[${i}] = ur[${i}];
    ur[${i}] = tmp;
    % endfor
</%pyfr:kernel>
