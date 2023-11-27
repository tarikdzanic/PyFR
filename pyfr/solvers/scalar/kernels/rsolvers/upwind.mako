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
        // Estimate the maximum wave speed 
        fpdtype_t lam = fmax(abs(laml), abs(lamr));
        nf[0] = 0.5*(fnl + fnr) + 0.5*lam*(ul[0] - ur[0]);
        % elif system == 'kpp':
        fpdtype_t laml = n[0]*cos(ul[0]) - n[1]*sin(ul[0]);
        fpdtype_t lamr = n[0]*cos(ur[0]) - n[1]*sin(ur[0]);

        fpdtype_t ap = max(0, max(laml, lamr));
        fpdtype_t am = min(0, min(laml, lamr));
        fpdtype_t da = ap - am;
        
        nf[0] = (ap*fnl - am*fnr)/da + (ap*am)*(ur[0] - ul[0])/da;

        % endif
    % endif
</%pyfr:macro>
