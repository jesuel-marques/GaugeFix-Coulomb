#ifndef SU2OPS_H
#define SU2OPS_H

void SU2_copy(const double* u, double* u_copy);

void SU2_set_to_null(double complex* u);
void SU2_set_to_identity(double complex* u);

void SU2_accumulate(const double* u, double* acc);

void SU2_subtraction(const double* u, const double* v, double* u_minus_v);

double SU2_determinant(const double* u);

void SU2_hermitean_conjugate(const double* u, double* u_dagger);

void SU2_multiplication_by_scalar(const double* u, const double alpha, double* alpha_times_u);

double SU2_inner_prod(const double* u, const double* v);

void SU2_outer_product(const double* u, const double* v, double* outer_product);
void SU2_product(const double* u, const double* v, double* uv);

void SU2_product_three(const double* u, const double* v, const double* w, double* uvw);

void SU2_product_four(const double* u, const double* v, const double* w, const double* x, double* uvwx);

void SU2_projection(double* a);

#endif