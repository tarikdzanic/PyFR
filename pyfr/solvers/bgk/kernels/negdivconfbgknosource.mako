# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfbgknosource' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

// Compute advection component df/dt = -u.div(f)
for (int i = 0; i < ${nvars}; i++) {
    tdivtconf[i] = -rcpdjac*tdivtconf[i];
}
</%pyfr:kernel>
