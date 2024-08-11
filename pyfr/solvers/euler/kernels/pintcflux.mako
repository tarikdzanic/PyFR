<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.euler.kernels.rotate'/>

<%pyfr:kernel name='pintcflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              nr='in fpdtype_t[${str(ndims)}]'
              rote='in fpdtype_t'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    // Normalize [nx,ny] vector only for rotations
    fpdtype_t mag_nl2 = sqrt(${pyfr.dot('nl[{i}]', i=2)});
    fpdtype_t norm_nl2[] = ${pyfr.array('(1 / mag_nl2)*nl[{i}]', i=2)};
    fpdtype_t negnorm_nr2[] = ${pyfr.array('-(1 / mag_nl2)*nr[{i}]', i=2)};

    ${pyfr.expand('rotate', 'ur', 'negnorm_nr2', 'norm_nl2')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'norm_nl', 'ficomm', 'rote')};

% for i in range(nvars):
    ul[${i}] =  mag_nl*(ficomm[${i}]);
    ur[${i}] = -mag_nl*(ficomm[${i}]);
% endfor

    // Transform RHS flux
    ${pyfr.expand('rotate', 'ur', 'norm_nl2', 'negnorm_nr2')};
</%pyfr:kernel>
