<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='negdivconfLO' ndim='1'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              rcpdjac='in fpdtype_t[${str(nupts)}]'
              ffpts='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'>


fpdtype_t nf[${nvars}], n[${ndims}], ul[${nvars}], ur[${nvars}];
fpdtype_t h2 = sqrt(1.0/rcpdjac[0]);

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = 0.0;
% endfor

n[0] = 1.0;
n[1] = 0.0;
% for i,j in pyfr.ndrange(p, p+1):
    // Create left/right states from u_i and u_i+1
    % for k in range(nvars):
        ul[${k}] = u[${(i)   + (p+1)*(j)}][${k}];
        ur[${k}] = u[${(i+1) + (p+1)*(j)}][${k}];
    % endfor
    // Compute f_i+1/2
    ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
    // Compute component of divF(u_i) and divF(u_+1)
    % for k in range(nvars):
        tdivtconf[${(i)   + (p+1)*(j)}][${k}] += ${ 1.0/(wts[i])}*nf[${k}]*h2;
        tdivtconf[${(i+1) + (p+1)*(j)}][${k}] += ${-1.0/(wts[i+1])}*nf[${k}]*h2;
    % endfor
% endfor

n[0] = 0.0;
n[1] = 1.0;
% for i,j in pyfr.ndrange(p+1, p):
    // Create left/right states from u_j and u_j+1
    % for k in range(nvars):
        ul[${k}] = u[${(i) + (p+1)*(j)  }][${k}];
        ur[${k}] = u[${(i) + (p+1)*(j+1)}][${k}];
    % endfor
    // Compute f_j+1/2
    ${pyfr.expand('rsolve', 'ul', 'ur', 'n', 'nf')};
    % for k in range(nvars):
        tdivtconf[${(i) + (p+1)*(j)  }][${k}] += ${ 1.0/(wts[j]  )}*nf[${k}]*h2;
        tdivtconf[${(i) + (p+1)*(j+1)}][${k}] += ${-1.0/(wts[j+1])}*nf[${k}]*h2;
    % endfor
% endfor
// Add face terms (assume symmetric weights)
% for fidx, k in pyfr.ndrange(p+1, nvars):
    tdivtconf[${(fidx) + (p+1)*(0)}][${k}]    += ffpts[${fidx}][${k}]*${1.0/wts[0]};
    tdivtconf[${(p)    + (p+1)*(fidx)}][${k}] += ffpts[${fidx + (p+1)}][${k}]*${1.0/wts[0]};
    tdivtconf[${(fidx) + (p+1)*(p)}][${k}]    += ffpts[${fidx + 2*(p+1)}][${k}]*${1.0/wts[0]};
    tdivtconf[${(0)    + (p+1)*(fidx)}][${k}] += ffpts[${fidx + 3*(p+1)}][${k}]*${1.0/wts[0]};
% endfor

% for i,j in pyfr.ndrange(nupts, nvars):
    tdivtconf[${i}][${j}] = -rcpdjac[0]*tdivtconf[${i}][${j}];
% endfor
## printf("%f %f %f %f \n", tdivtconf[0][0], tdivtconf[1][0], tdivtconf[2][0], tdivtconf[3][0]);
## printf("%f %f %f %f %f %f %f %f %f\n", tdivtconf[0][0], tdivtconf[1][0], tdivtconf[2][0], tdivtconf[3][0], tdivtconf[4][0], tdivtconf[5][0], tdivtconf[6][0], tdivtconf[7][0], tdivtconf[8][0]);

</%pyfr:kernel>
