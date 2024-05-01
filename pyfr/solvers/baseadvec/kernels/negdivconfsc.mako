<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:macro name='compute_indicator', params='u, alpha'>
    fpdtype_t ivar[${nupts}];

    // Compute indicator variable ivar = rho*P
% for i in range(nupts):
    % if ndims == 2:
    ivar[${i}] = ${c['gamma'] - 1}*(u[${i}][0]*u[${i}][${nvars-1}] - 0.5*(u[${i}][1]*u[${i}][1] + u[${i}][2]*u[${i}][2]));
    % elif ndims == 3:
    ivar[${i}] = ${c['gamma'] - 1}*(u[${i}][0]*u[${i}][${nvars-1}] - 0.5*(u[${i}][1]*u[${i}][1] + u[${i}][2]*u[${i}][2] + u[${i}][3]*u[${i}][3]));
    % endif
% endfor

    // Smoothness indicator
    fpdtype_t totEn = 0.0, pn1En = 1e-15, pn2En = 1e-15, tmp;
% for ivdm, b1mode, b2mode in zip(invvdm, ind1_modes, ind2_modes):
    tmp = ${' + '.join(f'{jx}*ivar[{j}]' for j, jx in enumerate(ivdm) if jx != 0)};

    totEn += tmp*tmp;
    % if b1mode:
    pn1En += tmp*tmp;
    % endif
    % if b2mode:
    pn2En += tmp*tmp;
    % endif
% endfor

    fpdtype_t En = fmax(pn1En/totEn, pn2En/totEn);

    alpha = 1.0/(1.0 + exp(-${s/Tn}*(En - ${Tn})));
    alpha = alpha < ${alpha_cutoff} ? 0 : alpha;
    alpha = alpha > ${1 - alpha_cutoff} ? 1 : alpha;
    alpha = fmin(alpha, ${alpha_max});
</%pyfr:macro>

<%pyfr:macro name='compute_tdivtconf_LO', params='u, rcpdjac, ffpts, pnx, pny, pnz, tdivtconf'>
    fpdtype_t nf[${nvars}], n[${ndims}], ul[${nvars}], ur[${nvars}], magn;

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = 0.0;
% endfor

% if ndims == 2:
    // i index changes first then j index
    <% uidx = lambda i,j : i + j*(p+1) %> 
    % for i,j in pyfr.ndrange(p, p+1):
        // Create left/right states from u_i and u_i+1
        % for k in range(nvars):
            ul[${k}] = u[${uidx(i,j)  }][${k}];
            ur[${k}] = u[${uidx(i+1,j)}][${k}];
        % endfor
        // Compute normal vector (averaged over left/right states for curved elements)
        % for k in range(ndims):
            n[${k}] = 0.5*(pnx[${uidx(i,j)}][${k}] + pnx[${uidx(i+1,j)}][${k}]);
        % endfor
        magn = sqrt(${pyfr.dot('n[{i}]', i=ndims)});
        % for k in range(ndims):
            n[${k}] /= magn;
        % endfor
        // Compute f_i+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        // Compute component of divF(u_i) and divF(u_i+1)
        % for k in range(nvars):
            tdivtconf[${uidx(i,j)  }][${k}] += ${ 1.0/(wts[i])  }*nf[${k}]*magn;
            tdivtconf[${uidx(i+1,j)}][${k}] += ${-1.0/(wts[i+1])}*nf[${k}]*magn;
        % endfor
    % endfor

    % for i,j in pyfr.ndrange(p+1, p):
        // Create left/right states from u_j and u_j+1
        % for k in range(nvars):
            ul[${k}] = u[${uidx(i,j)  }][${k}];
            ur[${k}] = u[${uidx(i,j+1)}][${k}];
        % endfor
        // Compute normal vector multiplied by face area (averaged over left/right states for curved elements)
        % for k in range(ndims):
            n[${k}] = 0.5*(pny[${uidx(i,j)}][${k}] + pny[${uidx(i,j+1)}][${k}]);
        % endfor
        magn = sqrt(${pyfr.dot('n[{i}]', i=ndims)});
        % for k in range(ndims):
            n[${k}] /= magn;
        % endfor
        // Compute f_j+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        % for k in range(nvars):
            tdivtconf[${uidx(i,j)  }][${k}] += ${ 1.0/(wts[j])  }*nf[${k}]*magn;
            tdivtconf[${uidx(i,j+1)}][${k}] += ${-1.0/(wts[j+1])}*nf[${k}]*magn;
        % endfor
    % endfor

    // Add face terms (assume symmetric weights: w[0] = w[-1])
    // Faces indices are bottom (-y), right (+x), top (+y), left (-x) for 2D
    <% fidx = lambda face, i : i + face*(p+1) %>
    % for idx, k in pyfr.ndrange(p+1, nvars):
        tdivtconf[${uidx(idx, 0)}][${k}] += ffpts[${fidx(0, idx)}][${k}]*${1.0/wts[0]}; // Bottom face
        tdivtconf[${uidx(p, idx)}][${k}] += ffpts[${fidx(1, idx)}][${k}]*${1.0/wts[0]}; // Right face
        tdivtconf[${uidx(idx, p)}][${k}] += ffpts[${fidx(2, idx)}][${k}]*${1.0/wts[0]}; // Top face
        tdivtconf[${uidx(0, idx)}][${k}] += ffpts[${fidx(3, idx)}][${k}]*${1.0/wts[0]}; // Left face
    % endfor
% elif ndims == 3:
    // i index changes first then j index then k
    <% uidx = lambda i,j,k : i + j*(p+1) + k*(p+1)**2 %>
    % for i,j,k in pyfr.ndrange(p, p+1, p+1):
        // Create left/right states from u_i and u_i+1
        % for l in range(nvars):
            ul[${l}] = u[${uidx(i,j,k)  }][${l}];
            ur[${l}] = u[${uidx(i+1,j,k)}][${l}];
        % endfor
        // Compute normal vector (averaged over left/right states for curved elements)
        % for l in range(ndims):
            n[${l}] = 0.5*(pnx[${uidx(i,j,k)}][${l}] + pnx[${uidx(i+1,j,k)}][${l}]);
        % endfor
        magn = sqrt(${pyfr.dot('n[{i}]', i=ndims)});
        % for l in range(ndims):
            n[${l}] /= magn;
        % endfor
        // Compute f_i+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        // Compute component of divF(u_i) and divF(u_i+1)
        % for l in range(nvars):
            tdivtconf[${uidx(i,j,k)  }][${l}] += ${ 1.0/(wts[i])  }*nf[${l}]*magn;
            tdivtconf[${uidx(i+1,j,k)}][${l}] += ${-1.0/(wts[i+1])}*nf[${l}]*magn;
        % endfor
    % endfor

    % for i,j,k in pyfr.ndrange(p+1, p, p+1):
        // Create left/right states from u_j and u_j+1
        % for l in range(nvars):
            ul[${l}] = u[${uidx(i,j,k)  }][${l}];
            ur[${l}] = u[${uidx(i,j+1,k)}][${l}];
        % endfor
        // Compute normal vector multiplied by face area (averaged over left/right states for curved elements)
        % for l in range(ndims):
            n[${l}] = 0.5*(pny[${uidx(i,j,k)}][${l}] + pny[${uidx(i,j+1,k)}][${l}]);
        % endfor
        magn = sqrt(${pyfr.dot('n[{i}]', i=ndims)});
        % for l in range(ndims):
            n[${l}] /= magn;
        % endfor
        // Compute f_j+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        % for l in range(nvars):
            tdivtconf[${uidx(i,j,k)  }][${l}] += ${ 1.0/(wts[j])  }*nf[${l}]*magn;
            tdivtconf[${uidx(i,j+1,k)}][${l}] += ${-1.0/(wts[j+1])}*nf[${l}]*magn;
        % endfor
    % endfor

    % for i,j,k in pyfr.ndrange(p+1, p+1, p):
        // Create left/right states from u_k and u_k+1
        % for l in range(nvars):
            ul[${l}] = u[${uidx(i,j,k)  }][${l}];
            ur[${l}] = u[${uidx(i,j,k+1)}][${l}];
        % endfor
        // Compute normal vector multiplied by face area (averaged over left/right states for curved elements)
        % for l in range(ndims):
            n[${l}] = 0.5*(pnz[${uidx(i,j,k)}][${l}] + pnz[${uidx(i,j,k+1)}][${l}]);
        % endfor
        magn = sqrt(${pyfr.dot('n[{i}]', i=ndims)});
        % for l in range(ndims):
            n[${l}] /= magn;
        % endfor
        // Compute f_k+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        % for l in range(nvars):
            tdivtconf[${uidx(i,j,k)  }][${l}] += ${ 1.0/(wts[k])  }*nf[${l}]*magn;
            tdivtconf[${uidx(i,j,k+1)}][${l}] += ${-1.0/(wts[k+1])}*nf[${l}]*magn;
        % endfor
    % endfor

    // Add face terms (assume symmetric weights: w[0] = w[-1])
    <% fidx = lambda face, i, j : i + j*(p+1) + face*(p+1)**2 %>
    // Faces indices are back (-z), bottom (-y), right (+x), top (+y), left (-x), front (+z) for 3D
    //       with i index changing first 
    % for iidx, jidx, k in pyfr.ndrange(p+1, p+1, nvars):
        tdivtconf[${uidx(iidx, jidx, 0)}][${k}] += ffpts[${fidx(0, iidx, jidx)}][${k}]*${1.0/wts[0]}; // Back face
        tdivtconf[${uidx(iidx, 0, jidx)}][${k}] += ffpts[${fidx(1, iidx, jidx)}][${k}]*${1.0/wts[0]}; // Bottom face
        tdivtconf[${uidx(p, iidx, jidx)}][${k}] += ffpts[${fidx(2, iidx, jidx)}][${k}]*${1.0/wts[0]}; // Right face
        tdivtconf[${uidx(iidx, p, jidx)}][${k}] += ffpts[${fidx(3, iidx, jidx)}][${k}]*${1.0/wts[0]}; // Top face
        tdivtconf[${uidx(0, iidx, jidx)}][${k}] += ffpts[${fidx(4, iidx, jidx)}][${k}]*${1.0/wts[0]}; // Left face
        tdivtconf[${uidx(iidx, jidx, p)}][${k}] += ffpts[${fidx(5, iidx, jidx)}][${k}]*${1.0/wts[0]}; // Front face
    % endfor
% endif

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = -rcpdjac[${i}]*tdivtconf[${i}][${j}];
% endfor
</%pyfr:macro>

              
<%pyfr:kernel name='negdivconfsc' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              ffpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              pnx='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              pny='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              pnz='in fpdtype_t[${str(nupts)}][${str(ndims)}]'>


% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = -rcpdjac[${i}]*tdivtconf[${i}][${j}];
% endfor

fpdtype_t alpha;
${pyfr.expand('compute_indicator', 'u', 'alpha')};

if (alpha > 0) {
    fpdtype_t tdivtconf_LO[${nupts}][${nvars}] = {{}};
    ${pyfr.expand('compute_tdivtconf_LO', 'u', 'rcpdjac', 'ffpts', 'pnx', 'pny', 'pnz', 'tdivtconf_LO')};

    % for i, j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = (1 - alpha)*tdivtconf[${i}][${j}] + alpha*tdivtconf_LO[${i}][${j}];
    % endfor
}
</%pyfr:kernel>
