<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
% for mod, name in src_macros:
    <%include file='${mod}'/>
% endfor

<%pyfr:kernel name='negdivconf' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
fpdtype_t src[${nvars}] = {};

% for mod, name in src_macros:
    ${pyfr.expand(name, 't', 'u', 'ploc', 'src')};
% endfor

% for i in range(nvars):
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + src[${i}];
% endfor

// Add rotational frame effects as source
// rhou_s = rho*omg*omg*x + 2*rhov*omg
// rhov_s = rho*omg*omg*y - 2*rhou*omg
// E_s    = -(mu*gamma/Pr)*(rote_xx + rote_yy) (assumes mu variation is negligible)
//        = -(mu*gamma/Pr)*(2*omg**2)
% if omg:
tdivtconf[1] += ${omg**2}*u[0]*ploc[0] + ${2*omg}*u[2];
tdivtconf[2] += ${omg**2}*u[0]*ploc[1] - ${2*omg}*u[1];
% if mu:
tdivtconf[${nvars-1}] += ${2*(mu*gamma/Pr)*omg**2};
% endif
% endif
</%pyfr:kernel>
