#include <tgmath.h>

#include <lapack.h>
#include <lapacke.h>

#include <SU3_ops.h>
#include <types.h>


static int eigensystem3x3(Mtrx3x3 * restrict a, 
                           Mtrx3x3 * restrict eigenvalues_mat, 
                           Mtrx3x3 * restrict eigenvectors) {
    double complex eigenvalues[Nc];
    
    Mtrx3x3 left_eigenvectors; /*not used*/

    lapack_int info;
    
    info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', Nc, a -> m, Nc, 
                        eigenvalues, 
                        &(left_eigenvectors).m, Nc, 
                             eigenvectors -> m, Nc);
    
    setNull3x3(eigenvalues_mat);
    MtrxIdx3 c;
    LOOP_3(c) {
        eigenvalues_mat -> m[ELEM_3X3(c, c)] = eigenvalues[c];
    }

    return info;
}


static void logDiagonalMtrx3x3(Mtrx3x3 * restrict diag_mtrx,  
                                  Mtrx3x3 * restrict log_of_diag_mtrx) {
    setNull3x3(log_of_diag_mtrx);    
    MtrxIdx3 a;
    LOOP_3(a) {
        log_of_diag_mtrx -> m[ELEM_3X3(a, a)] = log(diag_mtrx -> m[ELEM_3X3(a, a)]);
    }
}


static void powerDiagonalMtrx3x3(Mtrx3x3 * restrict diag_mtrx, Scalar power, 
                                    Mtrx3x3 * restrict diag_mtrx_to_power) {
    setNull3x3(diag_mtrx_to_power);    
    MtrxIdx3 a;
    LOOP_3(a) {
        diag_mtrx_to_power -> m[ELEM_3X3(a, a)] = pow(diag_mtrx -> m[ELEM_3X3(a, a)], 
                                                      power);
    }
}


int powerMtrx3x3(      Mtrx3x3 * restrict a, 
                     const Scalar power, 
                           Mtrx3x3 * restrict a_to_power) {
    Mtrx3x3 eigenvalues, eigenvalues_to_power;
    Mtrx3x3 eigenvectors, eigenvectors_inv;

    Mtrx3x3 a_copy; //because lapack overwrites it

    copy3x3(a, &a_copy);

    int status;
    status = eigensystem3x3(&a_copy, &eigenvalues, &eigenvectors);

    inverse3x3(&eigenvectors, &eigenvectors_inv);

    
    powerDiagonalMtrx3x3(&eigenvalues, power, &eigenvalues_to_power);

    prodThree3x3(&eigenvectors, &eigenvalues_to_power, &eigenvectors_inv, a_to_power);


    return status;
}


int logMtrx3x3(Mtrx3x3 * restrict a,  
               Mtrx3x3 * restrict log_of_a) {
    Mtrx3x3 eigenvalues, log_of_eigenvalues;
    Mtrx3x3 eigenvectors, eigenvectors_inv;

    Mtrx3x3 a_copy; //because lapack overwrites it

    copy3x3(a, &a_copy);

    int status;

    status = eigensystem3x3(&a_copy, &eigenvalues, &eigenvectors);

    inverse3x3(&eigenvectors, &eigenvectors_inv);

    
    logDiagonalMtrx3x3(&eigenvalues, &log_of_eigenvalues);

    prodThree3x3(&eigenvectors, &log_of_eigenvalues, &eigenvectors_inv, log_of_a);


    return status;
}