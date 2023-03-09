<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.scalar.kernels.flux'/>

<% tol = 1e-6 %>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    ${pyfr.expand('flux', 'ul', 'fl')};
    ${pyfr.expand('flux', 'ur', 'fr')};

    // Get left and right normal fluxes
    fpdtype_t fnl = ${pyfr.dot('n[{i}]', 'fl[{i}][0]', i=ndims)};
    fpdtype_t fnr = ${pyfr.dot('n[{i}]', 'fr[{i}][0]', i=ndims)};

    // Get the left and right wavespeeds
    fpdtype_t laml = abs(ul[0]) > ${tol} ? abs(fnl/ul[0]) : ${tol};
    fpdtype_t lamr = abs(ur[0]) > ${tol} ? abs(fnr/ur[0]) : ${tol};

    // Estimate the maximum wave speed / 2
    fpdtype_t lam = fmax(laml, lamr);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + 0.5*lam*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
