<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_entropy' params='u, d, p, e'>
    d = u[0];
    fpdtype_t rcpd = 1.0/d;
    fpdtype_t E = u[${nvars - 1}];

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*rcpd*(${pyfr.dot('u[{i}]', i=(1, ndims + 1))}));

    % if entropy_type == 'physical':
    // Compute specific physical entropy
    e = (d > 0 && p > 0) ? p*pow(rcpd, ${c['gamma']}) : ${fpdtype_max};
    % elif entropy_type == 'log':
    // Compute numerical entropy
    e = (d > 0 && p > 0) ? d*(log(p) - ${c['gamma']}*log(d)) : ${fpdtype_max};
    % endif
</%pyfr:macro>
