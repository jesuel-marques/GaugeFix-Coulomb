#include <SU2_ops.h>
#include <SU3_ops.h>
#include <fields.h>
#include <geometry.h>
#include <math.h>
#include <ranlux.h>
#include <stdio.h>
#include <tgmath.h>

#define POW2(a) ((a) * (a))
#define METROPOLIS_EPS 0.1

static inline void sumStaples(Mtrx3x3* restrict U,
                              const PosVec position,
                              const LorentzIdx mu,
                              Mtrx3x3* restrict staple_sum) {
    PosVec position_plusmu = getNeighbour(position, mu, FRONT);

    Mtrx3x3 u1a, u1b, u1c, staple1;
    Mtrx3x3 u2a, u2b, u2c, staple2;

    setNull3x3(staple_sum);

    LorentzIdx nu;
    LOOP_LORENTZ(nu) {
        if (mu != nu) {
            getLinkMatrix(U, position_plusmu, nu, FRONT, &u1a);
            getLinkMatrix(U, getNeighbour(position_plusmu, nu, FRONT), mu, REAR, &u1b);
            getLinkMatrix(U, getNeighbour(position, nu, FRONT), nu, REAR, &u1c);
            prodThree3x3(&u1a, &u1b, &u1c, &staple1);

            accumulate3x3(&staple1, staple_sum);

            getLinkMatrix(U, position_plusmu, nu, REAR, &u2a);
            getLinkMatrix(U, getNeighbour(position_plusmu, nu, REAR), mu, REAR, &u2b);
            getLinkMatrix(U, getNeighbour(position, nu, REAR), nu, FRONT, &u2c);
            prodThree3x3(&u2a, &u2b, &u2c, &staple2);

            accumulate3x3(&staple2, staple_sum);
        }
    }
}

static inline void HeatBathSubmatrix(Mtrx3x3* restrict u, Submtrx sub, Mtrx3x3* restrict W, double beta) {
    MtrxIdx3 a = sub == T ? 1 : 0;
    MtrxIdx3 b = sub == R ? 1 : 2;
    double lambda_squared;
    Mtrx2x2CK x;

    double modx, modf, mod_sqrd_f;

    Mtrx2x2CK w_dagger;
    Mtrx2x2CK update_SU2;

    w_dagger.m[0] = 0.5 * creal(W->m[ELEM_3X3(a, a)] + W->m[ELEM_3X3(b, b)]);
    w_dagger.m[1] = -0.5 * cimag(W->m[ELEM_3X3(a, b)] + W->m[ELEM_3X3(b, a)]);
    w_dagger.m[2] = -0.5 * creal(W->m[ELEM_3X3(a, b)] - W->m[ELEM_3X3(b, a)]);
    w_dagger.m[3] = -0.5 * cimag(W->m[ELEM_3X3(a, a)] - W->m[ELEM_3X3(b, b)]);

    double detw = projectSU2(&w_dagger);

    double e[4], f[3];

    do {
        ranlxd(e, 4);

        for (int i = 0; i < 4; i++) {
            e[i] = 1.0 - e[i];
        }

        lambda_squared = -((double)1.0 / (2.0 * (2.0 / 3.0) * sqrt(detw) * beta)) * (log(e[1]) + POW2(cos(2.0 * M_PI * e[2])) * log(e[3]));
    } while (POW2(e[0]) > 1.0 - lambda_squared);

    x.m[0] = 1.0 - 2.0 * lambda_squared;

    modx = sqrt(1.0 - POW2(x.m[0]));

    do {
        ranlxd(f, 3);

        for (int c = 0; c < 3; c++) {
            f[c] = 1.0 - 2.0 * f[c];
            mod_sqrd_f += POW2(f[c]);
        }

    } while (mod_sqrd_f > 1.0);

    modf = sqrt(mod_sqrd_f);

    for (int c = 1; c <= 3; c++) {
        x.m[c] = ((double)f[c - 1] / modf) * modx;
    }

    product2x2(&x, &w_dagger, &update_SU2);

    accumProdSU2_3x3(&update_SU2, u, a, b);
}

double HeatBathSU3(Mtrx3x3* restrict U, PosVec position, LorentzIdx mu, double beta) {
    // Atualização de elo utilizando algoritmo de banho térmico.

    Mtrx3x3 staple_sum;

    Mtrx3x3 W;

    Mtrx3x3* u = getLink(U, position, mu);

    sumStaples(U, position, mu, &staple_sum);
    prod3x3(u, &staple_sum, &W);
    double old_action_contribution = creal(trace3x3(&W));

    for (Submtrx sub = R; sub <= T; sub++) {
        HeatBathSubmatrix(u, sub, &W, beta);
        prod3x3(u, &staple_sum, &W);
    }

    projectSU3(u);

    double new_action_contribution = creal(trace3x3(&W));
    return -(beta / Nc) * (new_action_contribution - old_action_contribution);
}

static inline void MetropolisSubmatrix(Mtrx3x3* restrict update, Submtrx sub) {
    MtrxIdx3 a = sub == T ? 1 : 0;
    MtrxIdx3 b = sub == R ? 1 : 2;

    double r[4];
    ranlxd(r, 4);

    Mtrx2x2CK x;

    double norm_squared;

    for (int a = 1; a < 3; a++) {
        norm_squared += POW2(r[a]);
    }
    if (norm_squared != 0) {
        for (int a = 1; a < 3; a++) {
            r[a] /= sqrt(norm_squared);
        }
    } else {
        for (int a = 1; a < 3; a++) {
            r[a] = 1.0;
        }
    }

    x.m[0] = (r[0] > 0 ? 1.0 : -1.0) * sqrt(1.0 - POW2(METROPOLIS_EPS));
    for (int a = 1; a < 3; a++) {
        x.m[a] = METROPOLIS_EPS * r[a];
    }

    accumProdSU2_3x3(&x, update, a, b);
}

double MetropolisSU3(Mtrx3x3* restrict U, PosVec position, LorentzIdx mu, double beta) {
    Mtrx3x3 staple_sum;

    Mtrx3x3 W;

    Mtrx3x3* u = getLink(U, position, mu);

    Mtrx3x3 update, update_dagger;
    setIdentity3x3(&update);

    sumStaples(U, position, mu, &staple_sum);
    prod3x3(u, &staple_sum, &W);
    double old_action_contribution = creal(trace3x3(&W));

    for (Submtrx sub = R; sub <= T; sub++) {
        MetropolisSubmatrix(&update, sub);
    }

    double dagger_or_not;
    hermConj3x3(&update, &update_dagger);
    ranlxd(&dagger_or_not, 1);
    update = dagger_or_not > 0.0 ? update : update_dagger;

    accumLeftProd3x3(&update, &W);
    double would_be_action_contribution = creal(trace3x3(&W));
    double new_action_contribution;

    if ((would_be_action_contribution - old_action_contribution) < 0.0) {
        new_action_contribution = would_be_action_contribution;
        accumLeftProd3x3(&update, u);
    } else {
        double acceptQ;
        ranlxd(&acceptQ, 1);
        if (acceptQ <= exp(-(would_be_action_contribution - old_action_contribution))) {
            accumLeftProd3x3(&update, u);
            new_action_contribution = would_be_action_contribution;
        } else {
            new_action_contribution = old_action_contribution;
        }
    }

    projectSU3(u);

    return -(beta / Nc) * (new_action_contribution - old_action_contribution);
}