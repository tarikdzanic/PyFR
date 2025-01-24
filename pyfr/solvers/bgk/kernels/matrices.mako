# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_2x2inverse' params='M, Minv, invdet'>
invdet = M[0][0] * M[1][1] 
       - M[0][1] * M[1][0] ;
invdet = 1 / invdet;

Minv[0][0] = invdet *   ( M[1][1] );
Minv[0][1] = invdet * - ( M[0][1] );
Minv[1][0] = invdet * - ( M[1][0] );
Minv[1][1] = invdet *   ( M[0][0] );
</%pyfr:macro>

<%pyfr:macro name='compute_3x3inverse' params='M, Minv, invdet'>
invdet = M[0][0] * ( M[1][1] * M[2][2] - M[1][2] * M[2][1] ) 
       - M[0][1] * ( M[1][0] * M[2][2] - M[1][2] * M[2][0] ) 
       + M[0][2] * ( M[1][0] * M[2][1] - M[1][1] * M[2][0] ) ;
invdet = 1 / invdet;

Minv[0][0] = invdet *   ( M[1][1] * M[2][2] - M[1][2] * M[2][1] );
Minv[0][1] = invdet * - ( M[0][1] * M[2][2] - M[0][2] * M[2][1] );
Minv[0][2] = invdet *   ( M[0][1] * M[1][2] - M[0][2] * M[1][1] );
Minv[1][0] = invdet * - ( M[1][0] * M[2][2] - M[1][2] * M[2][0] );
Minv[1][1] = invdet *   ( M[0][0] * M[2][2] - M[0][2] * M[2][0] );
Minv[1][2] = invdet * - ( M[0][0] * M[1][2] - M[0][2] * M[1][0] );
Minv[2][0] = invdet *   ( M[1][0] * M[2][1] - M[1][1] * M[2][0] );
Minv[2][1] = invdet * - ( M[0][0] * M[2][1] - M[0][1] * M[2][0] );
Minv[2][2] = invdet *   ( M[0][0] * M[1][1] - M[0][1] * M[1][0] );
</%pyfr:macro>

<%pyfr:macro name='solve_linear_system' params='A, b, x, N'>
// N must be set outside of macro
fpdtype_t L[N][N] = {{0}};
fpdtype_t y[N], tmp;
int P[N];

// Set L to identity and initialize permutation vector
for (int i = 0; i < N; i++) {
    L[i][i] = 1.0;
    P[i] = i;
}

// Perform partial pivoting
for (int i = 0; i < N; i++) {
    // Find the row with the largest absolute value in column i
    int max_row = i;
    for (int j = i + 1; j < N; j++) {
        if (fabs(A[j][i]) > fabs(A[max_row][i])) {
            max_row = j;
        }
    }

    // Swap rows if needed
    if (max_row != i) {
        for (int k = 0; k < N; k++) {
            fpdtype_t temp = A[i][k];
            A[i][k] = A[max_row][k];
            A[max_row][k] = temp;
        }

        // Update the permutation vector
        int temp_idx = P[i];
        P[i] = P[max_row];
        P[max_row] = temp_idx;
    }
}

// Perform LU factorization (overwrite A with U)
for (int i = 0; i < N; i++) {
    for (int j = i+1; j < N; j++) {
        tmp = A[j][i]/A[i][i];
        L[j][i] = tmp;
        for (int k = 0; k < N; k++) {
            A[j][k] -= tmp*A[i][k];
        }
    }
}

// Get permuted b vector
fpdtype_t b_permuted[N];
for (int i = 0; i < N; i++) {
    b_permuted[i] = b[P[i]];
}

// Forward substitution
y[0] = b_permuted[0]/L[0][0];
for (int i = 1; i < N; i++) {
    tmp = 1.0/L[i][i];
    y[i] = tmp*b_permuted[i];
    for (int j = 0; j < i; j++) {
        y[i] -= tmp*L[i][j]*y[j];
    }
}

// Backward substitution
x[N-1] = y[N-1]/A[N-1][N-1];
for (int i = N-2; i > -1; i--) {
    tmp = 1/A[i][i];
    x[i] = tmp*y[i];
    for (int j = i+1; j < N; j++) {
        x[i] -= tmp*A[i][j]*x[j];
    }
}
</%pyfr:macro>