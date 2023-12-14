<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='correct' ndim='1'
              uHO='inout fpdtype_t[${str(nupts)}][${str(nvars)}]'
              uLO='in fpdtype_t[${str(nupts)}][${str(nvars)}]'>



fpdtype_t dminLO = ${fpdtype_max};
fpdtype_t dminHO = ${fpdtype_max};
fpdtype_t dmaxLO = ${-fpdtype_max};
fpdtype_t dmaxHO = ${-fpdtype_max};

% for i in range(nupts):
dminLO = fmin(dminLO, uLO[${i}][0]);
dminHO = fmin(dminHO, uHO[${i}][0]);
dmaxLO = fmax(dmaxLO, uLO[${i}][0]);
dmaxHO = fmax(dmaxHO, uHO[${i}][0]);
% endfor

if (dminHO < ${1 - salpha}*dminLO - ${1e-10} || dmaxHO > ${1 + salpha}*dmaxHO + ${1e-10}) {
    % for i,j in pyfr.ndrange(nupts, nvars):
    uHO[${i}][${j}] = uLO[${i}][${j}];
    % endfor
}
</%pyfr:kernel>
