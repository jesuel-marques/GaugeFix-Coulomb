#include <stdio.h>
#include <stdbool.h>

#include <flags.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <SU3_parameters.h>
#include <types.h>


extern short n_T;
extern short n_SPC;

inline bool position_valid(PosVec position){
    if (position.t < n_T   &&
        position.i < n_SPC &&
        position.j < n_SPC &&
        position.k < n_SPC    ){
        return true;
    }
    else{
        return false;
    }
}


inline bool position_mu_valid(PosVec position, 
                              LorentzIdx mu){
    if (position.t < n_T   &&
        position.i < n_SPC &&
        position.j < n_SPC &&
        position.k < n_SPC &&
       (mu == T_INDX || mu == X_INDX || 
        mu == Y_INDX || mu == Z_INDX )){
        return true;
    }
    return false;
}


PosVec assign_position(const PosIndex x, 
                       const PosIndex y, 
                       const PosIndex z, 
                       const PosIndex t) {
    /* assigns x, y, z and t to a position vector */
    PosVec position = {.i = x, .j = y, .k = z, .t = t};

    if(!position_valid(position)){
        fprintf(stderr, "Position {%d, %d, %d, %d} is invalid. "
                        "Returning origin instead\n.", x, y, z, t);
        return assign_position(0, 0, 0, 0);
    }

    return position;
}


void print_pos_vec(const PosVec pos) {
    /* prints a position to the screen*/
    printf("x: %hu y: %hu z: %hu t: %hu\n", pos.i, pos.j, pos.k, pos.t);
}


inline PosVec hop_pos_plus(const PosVec u, 
                           const LorentzIdx mu) {
    /* 	Calculates the position immediately forward
    	in the direction mu, taken into account the
    	periodic boundary conditions. */
    #ifdef CHECK_POSITION_BOUNDS
        if(!position_mu_valid(u, mu)){
            printf("Position: ");
            print_pos_vec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }
    #endif  //CHECK_POSITION_BOUNDS

    PosVec u_plus_muhat;

    unsigned short v;

    switch (mu) {
        case X_INDX:
            u_plus_muhat.i = ((v = u.i + 1) != n_SPC ? v : 0);
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = ((v = u.j + 1) != n_SPC ? v : 0);
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = ((v = u.k + 1) != n_SPC ? v : 0);
            u_plus_muhat.t = u.t;
            break;

        case T_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = ((v = u.t + 1) != n_T   ? v : 0);
            break;

    }

    return u_plus_muhat;
}


inline PosVec hop_pos_minus(const PosVec u, 
                            const LorentzIdx mu) {
    /* Calculates the position immediately behind
    in the direction mu, taken into account the
    periodic boundary conditions. */
    #ifdef CHECK_POSITION_BOUNDS
        if(!position_mu_valid(u, mu)){
            printf("Position: ");
            print_pos_vec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }
    #endif  //CHECK_POSITION_BOUNDS

    PosVec u_minus_muhat;

    short v ;

     switch (mu) {
        case X_INDX:
            u_minus_muhat.i = ((v = u.i - 1) != -1 ? v : n_SPC - 1);
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = ((v = u.j - 1) != -1 ? v : n_SPC - 1);
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = ((v = u.k - 1) != -1 ? v : n_SPC - 1);   
            u_minus_muhat.t = u.t;
            break;

        case T_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = ((v = u.t - 1) != -1 ?   v : n_T - 1);
            break;
    }

    return u_minus_muhat;
}