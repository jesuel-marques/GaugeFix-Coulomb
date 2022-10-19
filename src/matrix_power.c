#include <tgmath.h>

#include <lapack.h>
#include <lapacke.h>

#include <types.h>
#include <SU3_ops.h>

static int eigensystem_3x3(mtrx_3x3 * restrict a, mtrx_3x3 * restrict eigenvalues_mat, mtrx_3x3 * restrict eigenvectors){

    double complex eigenvalues[Nc];
    
    mtrx_3x3 left_eigenvectors; /*not used*/

    lapack_int info;
    
    info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', Nc, a -> m, Nc, eigenvalues, &(left_eigenvectors.m), Nc, eigenvectors -> m, Nc);
    
    set_null_3x3(eigenvalues_mat);
    for(int a = 0; a < Nc; a++){
        eigenvalues_mat -> m[ELM3x3(a, a)] = eigenvalues[a];
    }
    // // print_matrix_3x3(eigenvalues_mat, "eigenvalues interno", 16);
    // // print_matrix_3x3(eigenvalues_mat, "eigenvalues interno", 16);
    // copy_3x3(eigenvectors, eigenvectors_inv);
    // // print_matrix_3x3(eigenvalues_mat, "eigenvalues interno", 16);

    // // printf("addr eigenvalues:%p\n addr eigenvectors %p\n addr eigenvectors_inv %p", eigenvalues_mat, eigenvectors, eigenvectors_inv);
    // int ipiv[Nc] = {0, 0, 0};

    // info |= LAPACKE_zgetrf(LAPACK_ROW_MAJOR, Nc, Nc, eigenvectors_inv -> m, Nc, ipiv);
    // info |= LAPACKE_zgetri(LAPACK_ROW_MAJOR, Nc,     eigenvectors_inv -> m, Nc, ipiv);

    // // print_matrix_3x3(eigenvalues_mat, "eigenvalues interno", 16);
    // // print_matrix_3x3(eigenvectors, "eigenvectors interno", 16);
    // // print_matrix_3x3(eigenvectors_inv, "eigenvectors_inv interno", 16);

    return info;
}


static void log_diagonal_mtrx_3x3(mtrx_3x3 * restrict diag_mtrx,  
                                    mtrx_3x3 * restrict log_of_diag_mtrx){
    set_null_3x3(log_of_diag_mtrx);    
    for(int a = 0; a < Nc; a++)
        log_of_diag_mtrx -> m[ELM3x3(a, a)] = log(diag_mtrx -> m[ELM3x3(a, a)]);
    
}


static void power_diagonal_mtrx_3x3(mtrx_3x3 * restrict diag_mtrx, scalar power, 
                                    mtrx_3x3 * restrict diag_mtrx_to_power){
    set_null_3x3(diag_mtrx_to_power);    
    for(int a = 0; a < Nc; a++)
        diag_mtrx_to_power -> m[ELM3x3(a, a)] = pow(diag_mtrx -> m[ELM3x3(a, a)], power);
    
}

void matrix_power_3x3(mtrx_3x3 * restrict a, const scalar power, 
                            mtrx_3x3 * restrict a_to_power){

    mtrx_3x3 eigenvalues, eigenvalues_to_power;
    mtrx_3x3 eigenvectors, eigenvectors_inv;

    mtrx_3x3 a_copy; //because lapack overwrites it

    copy_3x3(a, &a_copy);

    int status;
    
    // print_matrix_3x3(a, "matrix", 18);

    status = eigensystem_3x3(&a_copy, &eigenvalues, &eigenvectors);

    inverse_3x3(&eigenvectors, &eigenvectors_inv);

    // print_matrix_3x3(&eigenvalues, "eigenvalues", 16);
    // print_matrix_3x3(&eigenvectors, "eigenvectors", 16);
    // print_matrix_3x3(&eigenvectors_inv, "eigenvectors_inv", 16);
    
    power_diagonal_mtrx_3x3(&eigenvalues, power, &eigenvalues_to_power);
    // print_matrix_3x3(&eigenvalues_to_power, "eigenvalues to power", 16);

    prod_three_3x3(&eigenvectors, &eigenvalues_to_power, &eigenvectors_inv, a_to_power);

    // print_matrix_3x3(a_to_power, "power of matrix", 18);

    return status;
}


int matrix_log_3x3(mtrx_3x3 * restrict a,  
                            mtrx_3x3 * restrict log_of_a){

    mtrx_3x3 eigenvalues, log_of_eigenvalues;
    mtrx_3x3 eigenvectors, eigenvectors_inv;

    mtrx_3x3 a_copy; //because lapack overwrites it

    copy_3x3(a, &a_copy);

    int status;
    
    // print_matrix_3x3(a, "matrix", 18);

    status = eigensystem_3x3(&a_copy, &eigenvalues, &eigenvectors);

    inverse_3x3(&eigenvectors, &eigenvectors_inv);

    // print_matrix_3x3(&eigenvalues, "eigenvalues", 16);
    // print_matrix_3x3(&eigenvectors, "eigenvectors", 16);
    // print_matrix_3x3(&eigenvectors_inv, "eigenvectors_inv", 16);
    
    log_diagonal_mtrx_3x3(&eigenvalues, &log_of_eigenvalues);
    // print_matrix_3x3(&eigenvalues_to_power, "eigenvalues to power", 16);

    prod_three_3x3(&eigenvectors, &log_of_eigenvalues, &eigenvectors_inv, log_of_a);

    // print_matrix_3x3(a_to_power, "power of matrix", 18);

    return status;
}