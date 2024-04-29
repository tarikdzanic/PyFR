<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

% for mod, name in src_macros:
    <%include file='${mod}'/>
% endfor

<%pyfr:macro name='calc_tdivtconf_LO', params='u, rcpdjac, ffpts, ploc, tdivtconf'>
    fpdtype_t nf[${nvars}], n[${ndims}], ul[${nvars}], ur[${nvars}];
    fpdtype_t h2 = sqrt(1.0/rcpdjac[0]);

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = 0.0;
% endfor

% if ndims == 2:
    <% uidx = lambda i,j : i + (p+1)*(j) %>
    <% fidx = lambda face, i : i + face*(p+1) %>
    n[0] = 1.0; n[1] = 0.0;
    % for i,j in pyfr.ndrange(p, p+1):
        // Create left/right states from u_i and u_i+1
        % for k in range(nvars):
            ul[${k}] = u[${uidx(i,j)  }][${k}];
            ur[${k}] = u[${uidx(i+1,j)}][${k}];
        % endfor
        // Compute f_i+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        // Compute component of divF(u_i) and divF(u_i+1)
        % for k in range(nvars):
            tdivtconf[${uidx(i,j)  }][${k}] += ${ 1.0/(wts[i])  }*nf[${k}]*h2;
            tdivtconf[${uidx(i+1,j)}][${k}] += ${-1.0/(wts[i+1])}*nf[${k}]*h2;
        % endfor
    % endfor

    n[0] = 0.0; n[1] = 1.0;
    % for i,j in pyfr.ndrange(p+1, p):
        // Create left/right states from u_j and u_j+1
        % for k in range(nvars):
            ul[${k}] = u[${uidx(i,j)  }][${k}];
            ur[${k}] = u[${uidx(i,j+1)}][${k}];
        % endfor
        // Compute f_j+1/2
        ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
        % for k in range(nvars):
            tdivtconf[${uidx(i,j)  }][${k}] += ${ 1.0/(wts[j])  }*nf[${k}]*h2;
            tdivtconf[${uidx(i,j+1)}][${k}] += ${-1.0/(wts[j+1])}*nf[${k}]*h2;
        % endfor
    % endfor

    // Add face terms (assume symmetric weights)
    % for idx, k in pyfr.ndrange(p+1, nvars):
        tdivtconf[${uidx(idx, 0)}][${k}] += ffpts[${fidx(0, idx)}][${k}]*${1.0/wts[0]}; // Bottom face
        tdivtconf[${uidx(p, idx)}][${k}] += ffpts[${fidx(1, idx)}][${k}]*${1.0/wts[0]}; // Right face
        tdivtconf[${uidx(idx, p)}][${k}] += ffpts[${fidx(2, idx)}][${k}]*${1.0/wts[0]}; // Top face
        tdivtconf[${uidx(0, idx)}][${k}] += ffpts[${fidx(3, idx)}][${k}]*${1.0/wts[0]}; // Left face
    % endfor
% endif

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = -rcpdjac[${i}]*tdivtconf[${i}][${j}]; /// CHANGE rcpdjac HEREEE
% endfor
</%pyfr:macro>

              
<%pyfr:kernel name='negdivconf' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              ploc='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              ffpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'>

fpdtype_t srci[${nvars}] = {}, ploci[${nvars}] = {}, ui[${nvars}] = {};
fpdtype_t src[${nupts}][${nvars}] = {{}};

fpdtype_t tdivtconf_HO[${nupts}][${nvars}] = {{}};
fpdtype_t tdivtconf_LO[${nupts}][${nvars}] = {{}};

% for i in range(nupts):
    % for j in range(ndims):
    ploci[${i}] = ploc[${i}][${j}];
    % endfor
    % for j in range(nvars):
    ui[${i}] = u[${i}][${j}];
    % endfor

    % for mod, name in src_macros:
    ${pyfr.expand(name, 't', 'ui', 'ploci', 'srci')};
    % endfor
    
    % for j in range(nvars):
    src[${i}][${j}] = srci[${j}];
    % endfor
% endfor

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf_HO[${i}][${j}] = -rcpdjac[${i}]*tdivtconf[${i}][${j}];
% endfor


${pyfr.expand('calc_tdivtconf_LO', 'u', 'rcpdjac', 'ffpts', 'ploc', 'tdivtconf')};

fpdtype_t alpha = 0.0;
% for i, j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = (1 - alpha)*tdivtconf_HO[${i}][${j}] + alpha*tdivtconf_LO[${i}][${j}] + src[${i}][${j}];
% endfor
</%pyfr:kernel>
