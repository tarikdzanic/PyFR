<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='elementmean' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ubar_int='out fpdtype_t[${str(nfaces)}][${str(nvars)}]'>

    fpdtype_t tmp;
% for i in range(nvars):
    tmp = ${' + '.join('{jx}*u[{j}][{i}]'.format(i=i, j=j, jx=jx)
                         for j, jx in enumerate(wts) if jx != 0)};
    % for j in range(nfaces):
    ubar_int[${j}][${i}] = tmp;
    % endfor
% endfor

</%pyfr:kernel>
