<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.scalar.kernels.flux'/>

<% tol = 1e-6 %>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf, ploc'>
    // Compute the left and right fluxes
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    ${pyfr.expand('flux', 'ul', 'fl', 'ploc')};
    ${pyfr.expand('flux', 'ur', 'fr', 'ploc')};

    // Get left and right normal fluxes
    fpdtype_t fnl = ${pyfr.dot('n[{i}]', 'fl[{i}][0]', i=ndims)};
    fpdtype_t fnr = ${pyfr.dot('n[{i}]', 'fr[{i}][0]', i=ndims)};

    % if system == 'advection':
    fpdtype_t cn = ${' + '.join(f'n[{j}]*{v[j]}' for j in range(ndims))};
    
    nf[0] = cn > 0.0 ? fnl : fnr;
    % else:
        % if system == 'burgers':
        fpdtype_t laml = ${' + '.join(f'n[{j}]*ul[{j}]' for j in range(ndims))};
        fpdtype_t lamr = ${' + '.join(f'n[{j}]*ur[{j}]' for j in range(ndims))};
        % elif system == 'kpp':
        fpdtype_t laml = n[0]*cos(ul[0]) + n[1]*sin(ul[0]);
        fpdtype_t lamr = n[0]*cos(ur[0]) + n[1]*sin(ur[0]);
        % endif
        // Estimate the maximum wave speed 
        fpdtype_t lam = fmax(abs(laml), abs(lamr));

        // Output
    % for i in range(nvars):
        nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                    for j in range(ndims))})
                + 0.5*lam*(ul[${i}] - ur[${i}]);
    % endfor
    % endif
</%pyfr:macro>
