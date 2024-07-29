# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='gmfrcptau' ndim='2'
              t='scalar fpdtype_t'
              f='inout fpdtype_t[${str(nvars)}]'
              g='in fpdtype_t[${str(nvars)}]'
              tau='in fpdtype_t'>

    for (int i = 0; i < ${nvars}; i++) {
        f[i] = (g[i] - f[i])/tau;
    }
</%pyfr:kernel>
