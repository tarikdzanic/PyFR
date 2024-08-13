# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='imexsolve' ndim='2'
              t='scalar fpdtype_t'
              z='scalar fpdtype_t'
              f='inout fpdtype_t[${str(nvars)}]'
              g='in fpdtype_t[${str(nvars)}]'
              tau='in fpdtype_t'>

    for (int i = 0; i < ${nvars}; i++) {
        f[i] = (f[i]*tau + z*g[i])/(tau + z);
    }
</%pyfr:kernel>
