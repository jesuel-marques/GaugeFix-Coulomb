#ifndef SU2OPS_H
#define SU2OPS_H

typedef unsigned short SU2_color_index;

void SU2_copy(const double* u, double* u_copy);

void SU2_set_to_null(double * u);
void SU2_set_to_identity(double * u);

void SU2_accumulate(const double* u, double* acc);

void SU2_subtraction(const double* u, const double* v, double* u_minus_v);

inline double SU2_determinant(const double* u);

void SU2_hermitean_conjugate(const double* u, double* u_dagger);

void SU2_multiplication_by_scalar(const double* u, const double alpha, double* alpha_times_u);

void SU2_product(const double* u, const double* v, double* uv);

void SU2_product_three(const double* u, const double* v, const double* w, double* uvw);

void SU2_product_four(const double* u, const double* v, const double* w, const double* x, double* uvwx);

inline void SU2_projection(double* u);

#endif