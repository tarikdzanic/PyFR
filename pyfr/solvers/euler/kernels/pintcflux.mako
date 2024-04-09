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

    fpdtype_t mag_nr = sqrt(${pyfr.dot('nr[{i}]', i=ndims)});
    fpdtype_t negnorm_nr[] = ${pyfr.array('-(1 / mag_nr)*nr[{i}]', i=ndims)};

    ${pyfr.expand('rotate', 'ur', 'negnorm_nr', 'norm_nl')};

    // Perform the Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'norm_nl', 'ficomm', 'rote')};

% for i in range(nvars):
    ul[${i}] =  mag_nl*(ficomm[${i}]);
    ur[${i}] = -mag_nl*(ficomm[${i}]);
% endfor

    // Transform RHS flux
    ${pyfr.expand('rotate', 'ur', 'norm_nl', 'negnorm_nr')};
</%pyfr:kernel>
