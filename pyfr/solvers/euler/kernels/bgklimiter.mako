<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% eps = 1e-12 %>
<%pyfr:kernel name='bgklimiter' ndim='1'
              u='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uf='in fpdtype_t[${str(nfpts)}][${str(nvars)}]'
              bounds='in fpdtype_t[${str(nvars*2)}]'
              m0='in broadcast fpdtype_t[${str(nfpts)}][${str(nupts)}]'>
           
      // Compute min/max of variables within element (and flux points)
      fpdtype_t umin[${nvars}];
      fpdtype_t umax[${nvars}];

% for i,j in pyfr.ndrange(nupts + nfpts, nvars):
      % if i == 0:
      umin[${j}] = u[${i}][${j}];
      umax[${j}] = u[${i}][${j}];
      % elif i < nupts:
      umin[${j}] = fmin(umin[${j}], u[${i}][${j}]);
      umax[${j}] = fmax(umax[${j}], u[${i}][${j}]);
      % else:
      umin[${j}] = fmin(umin[${j}], uf[${i - nupts}][${j}]);
      umax[${j}] = fmax(umax[${j}], uf[${i - nupts}][${j}]);
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
            % if i == 0:
            theta = fmin(theta, abs( (fmax(bounds[0], ${d_min}) - uavg[${i}])/(umin[${i}] - uavg[${i}]) ));
            % else:
            theta = fmin(theta, abs( (bounds[${i}]              - uavg[${i}])/(umin[${i}] - uavg[${i}]) ));
            % endif
      }
      if (abs(umax[${i}] - uavg[${i}]) > ${eps}) {
            theta = fmin(theta, abs( (bounds[${i+nvars}] - uavg[${i}])/(umax[${i}] - uavg[${i}]) ));
      }
      % endfor
      theta = fmax(0.0, theta);

      // Apply limiting
      % for i,j in pyfr.ndrange(nupts, nvars):
      u[${i}][${j}] = uavg[${j}] + theta*(u[${i}][${j}] - uavg[${j}]);
      % endfor

      // Recompute flux point values and check for minimum pressure bound
      fpdtype_t uf2[${nfpts}][${nvars}];
      for (int fidx = 0; fidx < ${nfpts}; fidx++) {
            for (int vidx = 0; vidx < ${nvars}; vidx++) {
                  uf2[fidx][vidx] = ${pyfr.dot('m0[fidx][{k}]', 'u[{k}][vidx]', k=nupts)};
            }
      }

      // Compute minimum pressure within element
      fpdtype_t pmin = ${fpdtype_max}, d, v[${ndims}], p, m2;
% for i in range(nupts + nfpts):
% if i < nupts:
      % if ndims == 2:
      p = ${c['gamma']-1}*(u[${i}][${nvars-1}] - (0.5/u[${i}][0])*(  u[${i}][1]*u[${i}][1]
                                                                   + u[${i}][2]*u[${i}][2]));
      % elif ndims == 3:
      p = ${c['gamma']-1}*(u[${i}][${nvars-1}] - (0.5/u[${i}][0])*(  u[${i}][1]*u[${i}][1]
                                                                   + u[${i}][2]*u[${i}][2]
                                                                   + u[${i}][3]*u[${i}][3]));
      % endif
% else:
      % if ndims == 2:
      p = ${c['gamma']-1}*(uf2[${i-nupts}][${nvars-1}] - (0.5/uf2[${i-nupts}][0])*(  uf2[${i-nupts}][1]*uf2[${i-nupts}][1]
                                                                                   + uf2[${i-nupts}][2]*uf2[${i-nupts}][2]));
      % elif ndims == 3:
      p = ${c['gamma']-1}*(uf2[${i-nupts}][${nvars-1}] - (0.5/uf2[${i-nupts}][0])*(  uf2[${i-nupts}][1]*uf2[${i-nupts}][1]
                                                                                   + uf2[${i-nupts}][2]*uf2[${i-nupts}][2]
                                                                                   + uf2[${i-nupts}][3]*uf2[${i-nupts}][3]));
      % endif
% endif
      pmin = fmin(pmin, p);
% endfor

      // Compute mean pressure
      % if ndims == 2:
      fpdtype_t pbar = ${c['gamma']-1}*(uavg[${nvars-1}] - (0.5/uavg[0])*(  uavg[1]*uavg[1]
                                                                          + uavg[2]*uavg[2]));
      % elif ndims == 3:
      fpdtype_t pbar = ${c['gamma']-1}*(uavg[${nvars-1}] - (0.5/uavg[0])*(  uavg[1]*uavg[1]
                                                                          + uavg[2]*uavg[2]
                                                                          + uavg[3]*uavg[3]));
      % endif

      // Apply limiting for pressure if necessary
      if (pmin < ${p_min} && abs(pmin - pbar) > ${eps}) {
            theta = fmin(1.0, abs( (${2*p_min} - pbar)/(pmin - pbar) ));
            % for i,j in pyfr.ndrange(nupts, nvars):
            u[${i}][${j}] = uavg[${j}] + theta*(u[${i}][${j}] - uavg[${j}]);
            % endfor
      }
</%pyfr:kernel>
