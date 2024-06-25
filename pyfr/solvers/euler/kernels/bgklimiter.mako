<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% eps = 1e-12 %>
<%pyfr:kernel name='bgklimiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              bounds='in fpdtype_t[${str(nvars*2)}]'>
      
      // Compute min/max of variables within element (and flux points)
      fpdtype_t umin[${nvars}] =  {${fpdtype_max}};
      fpdtype_t umax[${nvars}] = {-${fpdtype_max}};
% for i,j in pyfr.ndrange(nupts+nfpts, nvars):
      % if i < nupts:
      umin[${j}] = fmin(umin[${j}], u[${i}][${j}]);
      umax[${j}] = fmax(umax[${j}], u[${i}][${j}]);
      % else:
      umin[${j}] = fmin(umin[${j}], uf[${i}][${j}]);
      umax[${j}] = fmax(umax[${j}], uf[${i}][${j}]);
      % endif
% endfor

      // Compute mean modes
      fpdtype_t uavg[${nvars}];
      % for i in range(nvars):
      uavg[${i}] = ${' + '.join(f'{jx}*u[{j}][{i}]' 
                                    for j, jx in enumerate(meanwts) if jx != 0)};
      % endfor

      // Compute limiting factor
      fpdtype_t theta = 1.0;
      % for i in range(nvars):
      if (abs(umin[${i}] - uavg[${i}]) > ${eps}) {
            theta = fmin(theta, abs( (bounds[${i}  ] - uavg[${i}])/(umin[${i}] - uavg[${i}]) ));
      }
      if (abs(umax[${i}] - uavg[${i}]) > ${eps}) {
            theta = fmin(theta, abs( (bounds[${2*i}] - uavg[${i}])/(umax[${i}] - uavg[${i}]) ));
      }
      % endfor
      theta = fmax(0.0, theta);

      // Apply limiting
      % for i,j in pyfr.ndrange(nupts, nvars):
      u[${i}][${j}] = uavg[${j}] + theta*(u[${i}][${j}] - uavg[${j}]);
      % endfor


      // Check for minimum pressure bound
      fpdtype_t pmin = ${-fpdtype_max}, pbar = 0.0, d, v[${ndims}], p;
      % for i in range(nupts + nfpts):
      % if i < nupts:
            d = u[${i}][0];
            % for j in range(ndims):
            v[${j}] = u[${i}][${1+j}]/d;
            % endfor
            p = ${c['gamma']-1}*(u[${i}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{i}]', i=ndims)});
            pbar += ${meanwts[i]}*p;
      % else:
            d = uf[${i}][0];
            % for j in range(ndims):
            v[${j}] = uf[${i}][${1+j}]/d;
            % endfor
            p = ${c['gamma']-1}*(uf[${i}][${nvars-1}] - 0.5*d*${pyfr.dot('v[{i}]', i=ndims)});
      % endif
            pmin = fmin(pmin, p);
      % endfor

      // Apply limiting for pressure if necessary
      if (pmin < ${p_min} && (pmin - pbar) > ${eps}) {
            theta = fmin(1.0, abs( (${2*p_min} - pbar)/(pmin - pbar) ));
            % for i,j in pyfr.ndrange(nupts, nvars):
            u[${i}][${j}] = uavg[${j}] + theta*(u[${i}][${j}] - uavg[${j}]);
            % endfor
      }
</%pyfr:kernel>
