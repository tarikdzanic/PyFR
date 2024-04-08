<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

// Performs rotation from unit vector n to unit vector m for vector v
<%pyfr:macro name='rotate' params='v, n, m'>
    // Rotation matrix 
    // [ a   b
    //   c   d ]
    fpdtype_t a, b, c, d;
    a = n[0]*m[0] + n[1]*m[1];
    b = m[0]*n[1] - n[0]*m[1];
    c = n[0]*m[1] - m[0]*n[1];
    d = n[0]*m[0] + n[1]*m[1];
    // Offset by 1 (rotate only u/v velocity components)
    fpdtype_t v1 = v[1];
    fpdtype_t v2 = v[2];
    v[1] = a*v1 + b*v2;
    v[2] = c*v1 + d*v2;
</%pyfr:macro>
