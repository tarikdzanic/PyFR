<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.euler.kernels.rotate'/>

<%pyfr:kernel name='pintconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              urin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              urout='out view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'
              nr='in fpdtype_t[${str(ndims)}]'>

    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    fpdtype_t mag_nr = sqrt(${pyfr.dot('nr[{i}]', i=ndims)});
    fpdtype_t negnorm_nr[] = ${pyfr.array('-(1 / mag_nr)*nr[{i}]', i=ndims)};

    fpdtype_t tmpu[${nvars}];

    // ----- Assumes c['ldg-beta'] == 0.5 -----
    % for i in range(nvars):
    tmpu[${i}] = urin[${i}];
    % endfor
    ${pyfr.expand('rotate', 'tmpu', 'negnorm_nr', 'norm_nl')};
    % for i in range(nvars):
    ulout[${i}] = tmpu[${i}];
    % endfor
    // ----------------------------------------

</%pyfr:kernel>