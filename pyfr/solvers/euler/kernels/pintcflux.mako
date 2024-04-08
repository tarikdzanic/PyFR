<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='pintcflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              nr='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t mag_nr = sqrt(${pyfr.dot('nr[{i}]', i=ndims)});
    fpdtype_t norm_nr[] = ${pyfr.array('(1 / mag_nr)*nr[{i}]', i=ndims)};

    // Compute transformed solution
    fpdtype_t ul_t[${nvars}];
    fpdtype_t ur_t[${nvars}];

    // Density, energy, and w-velocity (if 3D) are not modified
    ul_t[0] = ul[0];
    ur_t[0] = ur[0];
    % if ndims == 3:
    ul_t[3] = ul[3];
    ur_t[3] = ur[3];
    % endif
    ul_t[${nvars-1}] = ul[${nvars-1}];
    ur_t[${nvars-1}] = ur[${nvars-1}];

    // Compute normal velocity with respect to outward facing normal
    ul_t[1] =  norm_nl[0]*ul[1] +  norm_nl[1]*ul[2];
    ur_t[1] = -norm_nr[0]*ur[1] + -norm_nr[1]*ur[2];
    // Compute tangential velocity with respect to outward facing normal
    ul_t[2] = -norm_nl[1]*ul[1] +  norm_nl[0]*ul[2];
    ur_t[2] =  norm_nr[1]*ur[1] + -norm_nr[0]*ur[2];
    // Create pseudo normal vector n = [1, 0, 0]
    fpdtype_t nn[${ndims}] = {0};
    nn[0] = 1.0;

    // Perform the Riemann solve
    fpdtype_t fn[${nvars}];
    ${pyfr.expand('rsolve', 'ul_t', 'ur_t', 'nn', 'fn')};

    // Scale, transform, and write out the common normal fluxes
% for i in range(nvars):
    % if i == 1:
    ul[${i}] = mag_nl*(norm_nl[0]*fn[1] - norm_nl[1]*fn[2]);
    ur[${i}] = mag_nr*(norm_nr[0]*fn[1] - norm_nr[1]*fn[2]);
    % elif i == 2:
    ul[${i}] = mag_nl*(norm_nl[1]*fn[1] + norm_nl[0]*fn[2]);
    ur[${i}] = mag_nr*(norm_nr[1]*fn[1] + norm_nr[0]*fn[2]);
    % else:
    ul[${i}] =  mag_nl*fn[${i}];
    ur[${i}] = -mag_nr*fn[${i}];
    % endif
% endfor
</%pyfr:kernel>
