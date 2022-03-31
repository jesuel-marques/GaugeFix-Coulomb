#ifndef SU3OPS_H
#define SU3OPS_H

void SU3_print_matrix(const double complex *u, const char *name);

void SU3_copy(const double complex *u, double complex *u_copy);

void SU3_set_to_null(double complex *u);

void SU3_set_to_identity(double complex *u);

void SU3_accumulate(const double complex *u, double complex *acc);

void SU3_subtraction(const double complex *u, const double complex *v, double complex *u_minus_v);

double complex SU3_trace(const double complex *u);

double complex SU3_determinant(const double complex *u);

void SU3_hermitean_conjugate(const double complex *u, double complex *u_dagger);

void SU3_multiplication_by_scalar(const double complex alpha, const double complex *u, double complex *alpha_times_u);
void SU3_substitution_multiplication_by_scalar(const double complex alpha, double complex *u);

void SU3_product(const double complex *u, const double complex *v, double complex *uv);

void SU3_product_three(const double complex *u, const double complex *v, const double complex *w, double complex *uvw);

void SU3_product_four(const double complex *u, const double complex *v, const double complex *w, const double complex *x, double complex *uvwx);

void SU3_accumulate_left_product(const double complex *g, double complex *acc_prod);

void SU3_accumulate_right_product(double complex *acc_prod, const double complex *g);

void SU3_projection(double complex *x);

void SU3_decompose_algebra(const double complex *a, double *a_components);

#endif