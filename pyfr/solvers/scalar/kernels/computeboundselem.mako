<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvec.kernels.computebounds'/>
<%include file='pyfr.solvers.scalar.kernels.elementbounds'/>

<%pyfr:kernel name='computeboundselem' ndim='1'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              xf='in broadcast fpdtype_t[${str(nfptsperface)}][${str(ndims-1)}]'
              bounds='inout fpdtype_t[2]'
              bounds_l='in fpdtype_t[${str(nfaces)}]'
              bounds_h='in fpdtype_t[${str(nfaces)}]'>

    for (int fidx = 0; fidx < ${nfaces}; fidx++) {
        bounds[0] = fmin(bounds[0], bounds_l[fidx]);
        bounds[1] = fmax(bounds[1], bounds_h[fidx]);
    }
</%pyfr:kernel>
