<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.bgk.kernels.util'/>

<%pyfr:kernel name='macrostate' ndim='2'
              f='in fpdtype_t[${str(nvars)}]'
              mvars='out fpdtype_t[${str(nmvars)}]'
              u='in broadcast fpdtype_t[${str(nuvars)}][${str(ndims)}]'
              M='in broadcast fpdtype_t[1][${str(nuvars)}]'>

    // Get macroscopic state
    ${pyfr.expand('compute_moments', 'f', 'u', 'M', 'mvars')};
</%pyfr:kernel>
