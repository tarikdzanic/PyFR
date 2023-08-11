# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='rsolve' params='fl, fr, n, nF, u'>
    fpdtype_t un = ${pyfr.dot('u[{i}]', 'n[{i}]', i=ndims)};    
    nF= un > 0.0 ? un*fl : un*fr;
</%pyfr:macro>