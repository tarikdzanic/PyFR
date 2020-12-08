# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='divf_LO_centered' ndim='1'
              u='in fpdtype_t[${str(nupts)}][${str(nvars)}]'
              plocu='in fpdtype_t[${str(nupts)}][${str(ndims)}]'
              usmats='in fpdtype_t[${str(nupts)}][${str(ndims*ndims)}]'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              fsmats='in fpdtype_t[${str(nfpts)}][${str(ndims*ndims)}]'
              divf='out fpdtype_t[${str(nupts)}][${str(nvars)}]'
              >

% for i,j in pyfr.ndrange(nupts, nvars):
	divf[${i}][${j}] = 0.0;
% endfor


fpdtype_t ftemp[${ndims}][${nvars}], ftemp2[${ndims}][${nvars}], fntemp[${nvars}], fntemp2[${nvars}, xl[${ndims}], xr[${ndims}];
fpdtype_t line_sol[${order+1}][${nvars}], line_flux[${order+2}][${ndims}][${nvars}], line_tflux[${order+2}][${ndims}][${nvars}];
fpdtype_t uk[${nvars}], n[${ndims}], t[${ndims}], tmp, p, v[${ndims}];

// Perform along xi direction
% for j in range(order+1):
	// Gather solution along constant eta line
	% for i in range(order+1):
		% for var in range(nvars):
			line_sol[${i}][${var}] = u[${i+ j*(order+1)}][${var}];
		% endfor
	% endfor

	% for i in range(order+1):
		// Set solution states
		% for var in range(nvars):
			uk[${var}] = line_sol[${i}][${var}];
		% endfor

        ${pyfr.expand('inviscid_flux', 'uk', 'ftemp', 'p', 'v')};
        % for dim, var in pyfr.ndrange(ndims, nvars):
            line_flux[${i+1}][${dim}][${var}] = ftemp[${dim}][${var}];
        % endfor

	% endfor

	// Definitely ndims*dim+k
	% for i in range(order+1):
		% for dim, var in pyfr.ndrange(ndims, nvars):
			line_tflux[${i}][${dim}][${var}] = ${' + '.join('(0.5*usmats[{0}][{2}] + 0.5*usmats[{1}][{2}])*line_flux[{5}][{3}][{4}]'
                                                 .format(i-1+j*(order+1),i+j*(order+1), ndims*dim+k, k, var, i) for k in range(ndims))};
		% endfor
	% endfor

	// Calculate df/dxi at solution points
	% for var in range(nvars):
		% for i in range(order+1):
			tmp =  ${' + '.join('{mx}*line_tflux[{m}][0][{var}]'.format(m=m, mx=mx, var=var)
                       for m, mx in enumerate(diffmatLO[i]) if mx != 0)};
			divf[${i+ j*(order+1)}][${var}] += tmp;
		% endfor
	% endfor 
% endfor

// Perform along eta direction
% for i in range(order+1):
	// Gather solution along constant xi line
	% for j in range(order+1):
		% for var in range(nvars):
			line_sol[${j}][${var}] = u[${i+ j*(order+1)}][${var}];
		% endfor
	% endfor

	// Flux split points

	// Take interior interface flux 
	% for var in range(nvars):
		ul[${var}] = uf[${i}][${var}];
		ur[${var}] = uf[${i + 2*(order+1)}][${var}];
	% endfor 
    ${pyfr.expand('inviscid_flux', 'ul', 'ftemp', 'p', 'v')};
    ${pyfr.expand('inviscid_flux', 'ur', 'ftemp2', 'p', 'v')};

	% for dim, var in pyfr.ndrange(ndims, nvars):
		line_flux[0][${dim}][${var}] = ftemp[${dim}][${var}];
		line_flux[${order+1}][${dim}][${var}] = ftemp2[${dim}][${var}];
	% endfor

	% for j in range(order):
		// Get normal direction between solution points
		% for dim in range(ndims):
			xl[${dim}] = plocu[${i+ (j+1)*(order+1)}][${dim}];
			xr[${dim}] = plocu[${i+ j*(order+1)}][${dim}];
		% endfor		
		${pyfr.expand('get_normal','xl', 'xr', 'n')}
		${pyfr.expand('get_tangent2d','xl', 'xr', 't')}

		// Set left and right solution states
		% for var in range(nvars):
			ul[${var}] = line_sol[${j}][${var}];
			ur[${var}] = line_sol[${j+1}][${var}];
		% endfor

        ${pyfr.expand('inviscid_flux', 'ul', 'ftemp', 'p', 'v')};
        ${pyfr.expand('inviscid_flux', 'ur', 'ftemp2', 'p', 'v')};
		% for dim, var in pyfr.ndrange(ndims, nvars):
            line_flux[${j+1}][${dim}][${var}] = 0.5*ftemp[${dim}][${var}] + 0.5*ftemp2[${dim}][${var}];
		% endfor
	% endfor

	// Transform flux to computational space
	% for dim, var in pyfr.ndrange(ndims, nvars):
		line_tflux[${0}][${dim}][${var}] = ${' + '.join('fsmats[{0}][{1}]*line_flux[{4}][{2}][{3}]'
                                                 .format(0*(order+1) + i, ndims*dim+k, k, var, 0) for k in range(ndims))};
		line_tflux[${order+1}][${dim}][${var}] = ${' + '.join('fsmats[{0}][{1}]*line_flux[{4}][{2}][{3}]'
                                                 .format(2*(order+1) + i, ndims*dim+k, k, var, order+1) for k in range(ndims))};
	% endfor

	% for j in range(1,order+1):
		% for dim, var in pyfr.ndrange(ndims, nvars):
			line_tflux[${j}][${dim}][${var}] = ${' + '.join('(0.5*usmats[{0}][{2}] + 0.5*usmats[{1}][{2}])*line_flux[{5}][{3}][{4}]'
                                                 .format(i+(j-1)*(order+1),i+j*(order+1), ndims*dim+k, k, var, j) for k in range(ndims))};
		% endfor
	% endfor

	// Calculate dg/deta at solution points
	% for var in range(nvars):
		% for j in range(order+1):
			tmp =  ${' + '.join('{mx}*line_tflux[{m}][1][{var}]'.format(m=m, mx=mx, var=var)
                       for m, mx in enumerate(diffmatLO[j]) if mx != 0)};
			divf[${i+ j*(order+1)}][${var}] += tmp; 
		% endfor
	% endfor

% endfor





</%pyfr:kernel>
