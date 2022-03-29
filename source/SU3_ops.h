#ifndef SU3OPS_H
#define SU3OPS_H

void SU3_print_matrix(double complex * u, char *name);

void SU3_copy(double complex * u, double complex * u_copy);

void SU3_set_to_null(double complex * u);

void SU3_set_to_identity(double complex * u );

void SU3_accumulate(double complex * u, double complex * acc);

void SU3_subtraction(double complex * u, double complex * v, double complex * u_minus_v);


double complex SU3_trace(double complex * u);

double complex SU3_determinant(double complex * u);

void SU3_hermitean_conjugate(double complex * u, double complex * u_dagger);


void SU3_multiplication_by_scalar(double complex alpha, double complex * u, double complex * alpha_times_u);
void SU3_substitution_multiplication_by_scalar(double complex alpha, double complex * u);

void SU3_product(double complex * u, double complex * v, double complex * uv);

void SU3_product_three(double complex * u, double complex * v, double complex * w, complex double * uvw);

void SU3_product_four(double complex * u, double complex * v, double complex * w, double complex * x, double complex * uvwx);

void SU3_accumulate_left_product(double complex * g, double complex * acc_prod);

void SU3_accumulate_right_product(double complex * acc_prod, double complex * g);

void SU3_projection(double complex * x);


void SU3_decompose_algebra(double complex * a, double * a_components);

#endif