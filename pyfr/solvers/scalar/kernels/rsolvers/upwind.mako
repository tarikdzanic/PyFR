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
        fpdtype_t fe[2];

        if (ul[0] < ur[0]) {
            fe[0] = min(fl[0][0], fr[0][0]);
            if (ur[0] - ul[0] > ${2*pi}) {
                fe[0] = -1;
            }
            else {
                fpdtype_t reml = fmod(ul[0], ${2*pi});
                fpdtype_t remr = fmod(ur[0], ${2*pi});
                if (remr > ${1.5*pi} && reml < ${1.5*pi}) {
                    fe[0] = -1;
                }
                else if (ur[0] - ul[0] > remr + ${0.5*pi}) {
                    fe[0] = -1;
                }
            }

        }
        else {
            fe[0] = max(fl[0][0], fr[0][0]);
            if (ul[0] - ur[0] > ${2*pi}) {
                fe[0] = 1;
            }
            else {
                fpdtype_t reml = fmod(ul[0], ${2*pi});
                fpdtype_t remr = fmod(ur[0], ${2*pi});
                if (remr > ${0.5*pi} && reml < ${0.5*pi}) {
                    fe[0] = 1;
                }
                else if (ul[0] - ur[0] > remr + ${1.5*pi}) {
                    fe[0] = 1;
                }
            }
        }


        if (ul[0] < ur[0]) {
            fe[1] = min(fl[1][0], fr[1][0]);
            if (ur[0] - ul[0] > ${2*pi}) {
                fe[1] = -1;
            }
            else {
                fpdtype_t reml = fmod(ul[0], ${2*pi});
                fpdtype_t remr = fmod(ur[0], ${2*pi});
                if (remr > ${1.0*pi} && reml < ${1.0*pi}) {
                    fe[1] = -1;
                }
                else if (ur[0] - ul[0] > remr + ${1.0*pi}) {
                    fe[1] = -1;
                }
            }

        }
        else {
            fe[1] = max(fl[1][0], fr[1][0]);
            if (ul[0] - ur[0] > ${2*pi}) {
                fe[1] = 1;
            }
            else {
                fpdtype_t reml = fmod(ul[0], ${2*pi});
                fpdtype_t remr = fmod(ur[0], ${2*pi});
                if (remr > ${0.5*pi} && reml < ${0.5*pi}) {
                    fe[1] = 1;
                }
                else if (reml < remr) {
                    fe[1] = 1;
                }
            }
        }

        fn[0] = n[0]*fe[0] + n[1]*fe[1];
    


        % endif
    % endif
</%pyfr:macro>
