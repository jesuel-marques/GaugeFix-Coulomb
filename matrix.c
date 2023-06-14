
#include <matrix.h>

#define DIGITS_MATRIX 18

/* Prints the 3x3 matrix to screen in Mathematica style. */
void printMatrix3x3(Mtrx u) {
    /*
     * Calls:
     * =====
     * printf,
     * creal, cimag.
     *
     * Macros:
     * ======
     * DIGITS_MATRIX,
     * Nc, ELEM_3X3, LOOP_3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:     3x3 matrix to be printed.
     *
     * Returns:
     * =======
     *
     */

    printf("{");
    MtrxIdx3 a, b;

    LOOP_MTRX(a, dim) {
        printf("{");
        LOOP_MTRX(b, dim) {
            printf("%.*lf + I*(%.*lf)", DIGITS_MATRIX, creal(u.*m[ELEM(a, b, dim)]),
                   DIGITS_MATRIX, cimag(u.*m[ELEM(a, b, dim)]));

            b != dim - 1 ? printf(", ") : 0;
        }
        a != dim - 1 ? printf("},\n") : 0;
    }
    printf("}}\n\n");
}

/* Copies a 3x3 matrix. */
void copy3x3(const Mtrx* restrict u, Mtrx* restrict u_copy) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix to be copied,
     * Mtrx3x3 * u_copy:    copied matrix.
     *
     * Returns:
     * =======
     *
     */

    memcpy(u_copy, u, sizeof(Mtrx3x3));
}

/* Sets a 3x3 matrix to be the 3x3 null matrix. */
inline void setNull3x3(Mtrx3x3* restrict u) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix to be set to null.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u->m[ELEM_3X3(a, b)] = 0.0;
    }
}

/* Sets a 3x3 matrix to be the 3x3 identity matrix. */
inline void setIdentity3x3(Mtrx3x3* restrict u) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix to be set to the 3x3 identity matrix.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u->m[ELEM_3X3(a, b)] = a != b ? 0.0 : 1.0;
    }
}

inline void setSU3Random(Mtrx3x3* u) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix to be set to the 3x3 identity matrix.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b;

    Mtrx2x2CK aux;
    setIdentity3x3(u);

    for (Submtrx sub = R; sub < T; sub++) {
        setSU2Random(&aux);

        a = sub == T ? 1 : 0;
        b = sub == R ? 1 : 2;
        accumProdSU2_3x3(&aux, u, a, b);
    }
}

/* Accumulates value of 3x3 matrix u in 3x3 matrix acc. */
inline void accumulate3x3(const Mtrx3x3* restrict u, Mtrx3x3* restrict acc) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix to be accumulated,
     * Mtrx3x3 * acc:       accumulated 3x3 matrix.
     *
     * Returns:
     * =======
     *
     */
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        acc->m[ELEM_3X3(a, b)] += u->m[ELEM_3X3(a, b)];
    }
}

/* Subtracts two 3x3 matrices and put value in another one. */
void subtraction3x3(Mtrx3x3* restrict u,
                    Mtrx3x3* restrict v,
                    Mtrx3x3* restrict u_minus_v, MtrxIndx dim) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix to be subtracted from,
     * Mtrx3x3 * v:         3x3 subtrahend matrix,
     * Mtrx3x3 * u_minus_v: difference matrix.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b;
    LOOP_MTRX(a, dim) {
        LOOP_MTRX(b) {
            u_minus_v->m[ELEM_3X3(a, b)] = u->m[ELEM_3X3(a, b)] - v->m[ELEM_3X3(a, b)];
        }
    }
}

/* Obtains the trace of a 3x3 matrix. */
Scalar trace(Mtrx* restrict u, MtrxIndx dim) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix whose trace is being calculated.
     *
     * Returns:
     * =======
     * Scalar containing the trace of the matrix,
     * which is a complex number in double precision.
     */

    double trace = 0.0;

    for (MtrxIdx a = 0; a < dim; a++) {
        trace += u->m[ELEM(a, a, dim)];
    }
    return trace;
}

/* Obtains the determinant of a 3x3 matrix. */
Scalar determinant3x3(Mtrx3x3* restrict u) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:     3x3 matrix whose determinant is being calculated.
     *
     * Returns:
     * =======
     * Scalar containing the determinant of the matrix,
     * which is a complex number in double precision.
     */

    //  The determinant is calculated in the usual manner  using Leibniz formula.
    return u->m[ELEM_3X3(0, 0)] * (u->m[ELEM_3X3(1, 1)] * u->m[ELEM_3X3(2, 2)] - u->m[ELEM_3X3(1, 2)] * u->m[ELEM_3X3(2, 1)]) + u->m[ELEM_3X3(0, 1)] * (u->m[ELEM_3X3(1, 2)] * u->m[ELEM_3X3(2, 0)] - u->m[ELEM_3X3(1, 0)] * u->m[ELEM_3X3(2, 2)]) + u->m[ELEM_3X3(0, 2)] * (u->m[ELEM_3X3(1, 0)] * u->m[ELEM_3X3(2, 1)] - u->m[ELEM_3X3(1, 1)] * u->m[ELEM_3X3(2, 0)]);
}

/* Calculates the hermitean conjugate of a 3x3 matrix. */
void hermConj(Mtrx* restrict u, Mtrx* restrict u_dagger, MtrxIndx dim) {
    /*
     * Calls:
     * =====
     * conj.
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:         3x3 matrix whose hermitean conjugate is being calculated,
     * Mtrx3x3 * u_dagger:  hermitean conjugated u.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx a, b;
    LOOP_MTRX(a) {
        LOOP_MTRX(b) {
            u_dagger->m[ELEM(b, a, dim)] = conj(u->m[ELEM(a, b, dim)]);
        }
    }
}

/* Calculates multiplication of 3x3 matrix u by Scalar alpha, which is a complex
   number in double precision. */
void multByScalar3x3(const Scalar alpha,
                     Mtrx3x3* restrict u,
                     Mtrx3x3* restrict alpha_times_u,
                     MtrxIndx dim) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Scalar alpha:            scalar to multiply the matrix,
     * Mtrx3x3 * u:             3x3 matrix to be multiplied by alpha,
     * Mtrx3x3 * alpha_times_u: matrix resulting of multiplication of alpha by u.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx a, b;
    LOOP_MTRX(a) {
        LOOP_MTRX(b) {
            alpha_times_u->m[ELEM(a, b, dim)] = alpha * u->m[ELEM(a, b, dim)];
        }
    }
}

/* Calculates multiplication of 3x3 matrix u by Scalar alpha, which is a complex
   number in double precision. Result is returned in same matrix. */
inline void substMultScalar3x3(const Scalar alpha,
                               Mtrx3x3* restrict u) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Scalar alpha:        scalar to multiply the matrix,
     * Mtrx3x3 * u:         3x3 matrix to be multiplied by alpha
     *                      and in which the result will be put.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u->m[ELEM_3X3(a, b)] *= alpha;
    }
}

/* Calculates product of two 3x3 matrices. */
inline void prod3x3(const Mtrx3x3* restrict u,
                    const Mtrx3x3* restrict v,
                    Mtrx3x3* restrict uv) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:     first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:     second 3x3 matrix to be taken in the product,
     * Mtrx3x3 * uv:    result of product.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b, c;

    LOOP_3X3(a, b) {
        uv->m[ELEM_3X3(a, b)] = 0.0;
        LOOP_3(c) {
            uv->m[ELEM_3X3(a, b)] += (u->m[ELEM_3X3(a, c)]) * (v->m[ELEM_3X3(c, b)]);
        }
    }
}

/* Calculates product of three 3x3 matrices. */
void prodThree3x3(const Mtrx3x3* restrict u,
                  const Mtrx3x3* restrict v,
                  const Mtrx3x3* restrict w,
                  Mtrx3x3* restrict uvw) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:     first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:     second 3x3 matrix to be taken in the product,
     * Mtrx3x3 * w:     third 3x3 matrix to be taken in the product,
     * Mtrx3x3 * uvw:   result of product.
     *
     * Returns:
     * =======
     *
     */

    Mtrx3x3 uv;
    prod3x3(u, v, &uv);
    prod3x3(&uv, w, uvw);
}

/* Calculates product of four 3x3 matrices. */
void prodFour3x3(const Mtrx3x3* restrict u,
                 const Mtrx3x3* restrict v,
                 const Mtrx3x3* restrict w,
                 const Mtrx3x3* restrict x,
                 Mtrx3x3* restrict uvwx) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:     first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:     second 3x3 matrix to be taken in the product,
     * Mtrx3x3 * w:     third 3x3 matrix to be taken in the product,
     * Mtrx3x3 * x:     fourth 3x3 matrix to be taken in the product,
     * Mtrx3x3 * uvwx:  result of product.
     *
     * Returns:
     * =======
     *
     */

    Mtrx3x3 aux1, aux2;

    prod3x3(u, v, &aux1);
    prod3x3(w, x, &aux2);
    prod3x3(&aux1, &aux2, uvwx);
}

void prod3x3_generic(Mtrx3x3** restrict u, short number_of_matrices,
                     Mtrx3x3* restrict result) {
    setIdentity3x3(result);
    for (int i = 0; i < number_of_matrices; i++) {
        accumRightProd3x3(result, *(u + i));
    }
}
inline void accumLeftProd3x3(const Mtrx3x3* restrict g,
                             Mtrx3x3* restrict acc_prod) {
    /*
     * Calls:
     * =====
     * copy3x3, prod3x3.
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * g:         3x3 matrix to be left product accumulated in the acc_prod
     * Mtrx3x3 * acc_prod:  3x3 matrix which holds the accumulated product
     *
     * Returns:
     * =======
     *
     */

    Mtrx3x3 aux_prod;

    copy3x3(acc_prod, &aux_prod);
    prod3x3(g, &aux_prod, acc_prod);
}

/* Calculates matrix product between acc_prod and g and accumulates result in acc_prod.
 */
inline void accumRightProd3x3(Mtrx3x3* restrict acc_prod,
                              const Mtrx3x3* restrict g) {
    /*
     * Calls:
     * =====
     * copy3x3, prod3x3.
     *
     * Macros:
     * ======
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * acc_prod:  3x3 matrix which holds the accumulated product
     * Mtrx3x3 * g:         3x3 matrix to be right product accumulated in the acc_prod
     *
     * Returns:
     * =======
     *
     */

    Mtrx3x3 aux_prod;

    copy3x3(acc_prod, &aux_prod);
    prod3x3(&aux_prod, g, acc_prod);
}

/* Calculates matrix product between v, u and w dagger. */
inline void prod_vuwdagger3x3(const Mtrx3x3* restrict v,
                              const Mtrx3x3* restrict u,
                              const Mtrx3x3* restrict w,
                              Mtrx3x3* restrict vuwdagger) {
    /*
     * Calls:
     * =====
     * conj
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * u:             first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:             middle 3x3 matrix to be taken in the product,
     * Mtrx3x3 * w:             3x3 matrix to be taken the hermitean conjugate and then,
     *                          taken in the product,
     * Mtrx3x3 * vuwdagger:     result of product.
     *
     * Returns:
     * =======
     *
     */

    MtrxIdx3 a, b, c, d;
    LOOP_3X3(a, b) {
        vuwdagger->m[ELEM_3X3(a, b)] = 0.0;
        LOOP_3(c) {
            LOOP_3(d) {
                vuwdagger->m[ELEM_3X3(a, b)] += (v->m[ELEM_3X3(a, c)]) * (u->m[ELEM_3X3(c, d)]) * conj(w->m[ELEM_3X3(b, d)]);
                //  wdagger_db is conj(w_bd)
            }
        }
    }
}

/* Calculates A^omega according to the formula (which works for A close to identity)
   A^omega = Proj_SU3(I + omega * (A-I) + ...). */
inline short power3x3Binomial(const Mtrx3x3* restrict A,
                              const double omega,
                              Mtrx3x3* restrict A_to_omega) {
    /*
     * Calls:
     * =====
     * projectSU3
     *
     * Macros:
     * ======
     * ELEM_3X3, LOOP_3X3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * A:             3x3 matrix whose power is being taken
     * double omega:            the power that is being taken
     * Mtrx3x3 * A_to_omega:    result of the exponentiation
     *
     * Returns:
     * =======
     *
     */

    //  A^omega = Proj_SU3(I + omega * (A-I) + ...)
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        A_to_omega->m[ELEM_3X3(a, b)] = a == b ? 1.0 + omega * (A->m[ELEM_3X3(a, b)] - 1.0) : omega * (A->m[ELEM_3X3(a, b)]);
    }

    return projectSU3(A_to_omega);
}

/* Takes a color vector and normalizes it. */
inline static int normalize3Vec(Vec3* restrict v) {
    /*
     * Calls:
     * =====
     * fabs, sqrt.
     *
     * Macros:
     * ======
     * POW2,
     * LOOP_3.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * v:     color vector to be normalized.
     *
     * Returns:
     * =======
     * 0 if successful, -1 otherwise. Normalization will not be successful if vector was
     * the null vector.
     */

    double r = 0.0;

    MtrxIdx3 a;

    LOOP_3(a) {
        //  Calculates the sum of the absolute
        //  value squared of each component
        r += (fabs(v->m[a])) * (fabs(v->m[a]));
    }
    if (r != 0) {
        //  If vector is not null, then normalize it
        r = 1.0 / sqrt(r);

        LOOP_3(a) {
            v->m[a] *= r;
        }
        return 0;
    } else {
        return -1;
    }
}
