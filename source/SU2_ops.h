#ifndef SU2OPS_H
#define SU2OPS_H

void SU2_copy(double* u, double* u_copy);

void SU2_set_to_null(double complex* u);
void SU2_set_to_identity(double complex* u);

void SU2_accumulate(double* u, double* acc);

void SU2_subtraction(double* u, double* v, double* u_minus_v);

double SU2_determinant(double* u);

void SU2_hermitean_conjugate(double* u, double* u_dagger);

void SU2_multiplication_by_scalar(double* u, double alpha, double* alpha_times_u);

double SU2_inner_prod(double* u, double* v);

void SU2_outer_product(double* u, double* v, double* outer_product);
void SU2_product(double* u, double* v, double* uv);

void SU2_product_three(double* u, double* v, double* w, double* uvw);

void SU2_product_four(double* u, double* v, double* w, double* x, double* uvwx);

void SU2_projection(double* a);

#endif