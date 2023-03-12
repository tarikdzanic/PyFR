<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='flux' params='u, f, ploc'>
    // Compute the flux
% if system == 'advection':
    % for i in range(ndims):
    f[${i}][0] = ${v[i]}*u[0];
    % endfor
% elif system == 'burgers':
    % for i in range(ndims):
    f[${i}][0] = ${v[i]}*u[0]*u[0];
    % endfor
% elif system == 'kpp':
    f[0][0] = sin(u[0]);
    f[1][0] = cos(u[0]);
% endif
</%pyfr:macro>
