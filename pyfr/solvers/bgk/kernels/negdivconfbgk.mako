# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='negdivconfbgk' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t'
              coll='in fpdtype_t'
              rcpdjac='in fpdtype_t'>

// Compute df/dt = -u.div(f) + (g-f)/tau
tdivtconf = -rcpdjac*tdivtconf + coll;
</%pyfr:kernel>
