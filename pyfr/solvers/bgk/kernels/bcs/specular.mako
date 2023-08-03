# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ntol = 1e-6 %>
<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr' externs='ploc, t'>

% if ndims == 2:
// If +/- X normal:
if (abs(abs(nl[0]) - 1.0) < ${ntol} && abs(nl[1]) < ${ntol}) {
	% for i in range(nvars):
	fr[${i}] = fl[${Xidxs[i]}];
	% endfor
}
// If +/- Y normal:
else if (abs(abs(nl[1]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol}) {
	% for i in range(nvars):
	fr[${i}] = fl[${Yidxs[i]}];
	% endfor
}
% elif ndims == 3:
// If +/- X normal:
if (abs(abs(nl[0]) - 1.0) < ${ntol} && abs(nl[1]) < ${ntol} && abs(nl[2]) < ${ntol}) {
	% for i in range(nvars):
	fr[${i}] = fl[${Xidxs[i]}];
	% endfor
}
// If +/- Y normal:
else if (abs(abs(nl[1]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol} && abs(nl[2]) < ${ntol}) {
	% for i in range(nvars):
	fr[${i}] = fl[${Yidxs[i]}];
	% endfor
}
// If +/- Z normal:
else if (abs(abs(nl[2]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol} && abs(nl[1]) < ${ntol}) {
	% for i in range(nvars):
	fr[${i}] = fl[${Zidxs[i]}];
	% endfor
}
% endif
</%pyfr:macro>
