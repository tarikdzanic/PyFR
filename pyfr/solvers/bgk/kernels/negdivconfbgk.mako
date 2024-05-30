# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfbgk' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              coll='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>

// Compute df/dt = -u.div(f) + (g-f)/tau
for (int i = 0; i < ${nvars}; i++) {
    tdivtconf[i] = -rcpdjac*tdivtconf[i] + coll[i];
}
</%pyfr:kernel>
