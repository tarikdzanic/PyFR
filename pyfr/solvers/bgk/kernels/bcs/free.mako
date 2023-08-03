# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr' externs='ploc, t'>
    for (int i = 0; i < ${nvars}; i++) {
        fr[i] = fl[i];
    }
</%pyfr:macro>
