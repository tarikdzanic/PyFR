<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.computebounds'/>
<%include file='pyfr.solvers.scalar.kernels.elementbounds'/>

<%pyfr:kernel name='computeboundsface' ndim='1'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              xf='in broadcast fpdtype_t[${str(nfptsperface)}][${str(ndims-1)}]'
              bounds='inout fpdtype_t[2]'>

    fpdtype_t ufi[${nfptsperface}][${nvars}], alpha;
    % for i in range(nfaces):
    // Get face nodal solution
    % for j,k in pyfr.ndrange(nfptsperface, nvars):
    ufi[${j}][${k}] = uf[${i*nfptsperface + j}][${k}];
    % endfor

    ${pyfr.expand('optimize_face_l', 'ufi' ,'xf', 'alpha')};
    bounds[0] = fmin(bounds[0], alpha);

    ${pyfr.expand('optimize_face_h', 'ufi' ,'xf', 'alpha')};
    bounds[1] = fmax(bounds[1], -alpha);

    ${pyfr.expand('optimize_face_e', 'ufi' ,'xf', 'alpha')};
    bounds[2] = fmin(bounds[2], alpha);
    % endfor
</%pyfr:kernel>
