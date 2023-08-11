# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<% ntol = 1e-6 %>
<%pyfr:macro name='bc_rsolve_state' params='fl, nl, fr, u, M' externs='ploc, t'>

% if ndims == 2:
// If +/- X normal:
if (abs(abs(nl[0]) - 1.0) < ${ntol} && abs(nl[1]) < ${ntol}) {
	% for i in range(nuvars):
	fr[${i}] = fl[${Xidxs[i]}];
	% endfor
	% if delta:
	% for i in range(nuvars):
	fr[${i + nuvars}] = fl[${Xidxs[i] + nuvars}];
	% endfor
	% endif
}
// If +/- Y normal:
else if (abs(abs(nl[1]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol}) {
	% for i in range(nuvars):
	fr[${i}] = fl[${Yidxs[i]}];
	% endfor
	% if delta:
	% for i in range(nuvars):
	fr[${i + nuvars}] = fl[${Yidxs[i] + nuvars}];
	% endfor
	% endif
}
% elif ndims == 3:
// If +/- X normal:
if (abs(abs(nl[0]) - 1.0) < ${ntol} && abs(nl[1]) < ${ntol} && abs(nl[2]) < ${ntol}) {
	% for i in range(nuvars):
	fr[${i}] = fl[${Xidxs[i]}];
	% endfor
	% if delta:
	% for i in range(nuvars):
	fr[${i + nuvars}] = fl[${Xidxs[i] + nuvars}];
	% endfor
	% endif
}
// If +/- Y normal:
else if (abs(abs(nl[1]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol} && abs(nl[2]) < ${ntol}) {
	% for i in range(nuvars):
	fr[${i}] = fl[${Yidxs[i]}];
	% endfor
	% if delta:
	% for i in range(nuvars):
	fr[${i + nuvars}] = fl[${Yidxs[i] + nuvars}];
	% endfor
	% endif
}
// If +/- Z normal:
else if (abs(abs(nl[2]) - 1.0) < ${ntol} && abs(nl[0]) < ${ntol} && abs(nl[1]) < ${ntol}) {
	% for i in range(nuvars):
	fr[${i}] = fl[${Zidxs[i]}];
	% endfor
	% if delta:
	% for i in range(nuvars):
	fr[${i + nuvars}] = fl[${Zidxs[i] + nuvars}];
	% endfor
	% endif
}
% endif
</%pyfr:macro>
