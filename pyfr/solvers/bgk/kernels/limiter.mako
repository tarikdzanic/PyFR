# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%pyfr:kernel name='limiter' ndim='1'
              f='inout fpdtype_t[${str(nupts)}]'>
<% eps = 1e-12%>

fpdtype_t avg_f = ${' + '.join('{jx}*f[{j}]'.format(j=j, jx=jx)
                        for j, jx in enumerate(wts) if jx != 0)};

fpdtype_t min_f = f[0];
% for j in range(1, nupts):
min_f = fmin(min_f, f[${j}]);
% endfor

if (min_f < 0) {
    fpdtype_t beta = abs(avg_f/max(avg_f - min_f, ${eps}));    
    beta = fmax(0.0, fmin(beta, 1.0));       
    % for j in range(nupts):
    f[${j}] = beta*(f[${j}] - avg_f) + avg_f;
    % endfor
}
</%pyfr:kernel>
