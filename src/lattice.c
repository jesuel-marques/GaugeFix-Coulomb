#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <complex.h>

#include <SU3_parameters.h>              //	Simulation parameters
#include <SU3_gaugefixing_parameters.h>  //	Gauge-fixing specific parameters

#include <SU2_ops.h>                        //	SU(2) operations
#include <SU3_ops.h>                        //	SU(3) operations
                                            //	positions and links on the lattice.

#include <lattice.h>  //	Initialization functions and calculations of

inline bool position_valid(pos_vec position){
    if (position.t >= 0 && position.t < N_T   &&
        position.i >= 0 && position.i < N_SPC &&
        position.j >= 0 && position.j < N_SPC &&
        position.k >= 0 && position.k < N_SPC    ){
        return true;
    }
    else{
        printf("Position invalid: ");
        print_pos_vec(position);
        printf("\n");
        return false;
    }
}

inline bool position_mu_valid(pos_vec position, lorentz_idx mu){
    if (position.t >= 0 && position.t < N_T   &&
        position.i >= 0 && position.i < N_SPC &&
        position.j >= 0 && position.j < N_SPC &&
        position.k >= 0 && position.k < N_SPC &&
       (mu == T_INDX || mu == X_INDX || 
        mu == Y_INDX || mu == Z_INDX )){

        return true;

    }

    return false;

}

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t) {
    //	assigns x, y, z and t to a position vector
    pos_vec position;

    position.i = x;
    position.j = y;
    position.k = z;
    position.t = t;

    if(!position_valid(position)){
        fprintf(stderr, "Position {%d, %d, %d, %d} is invalid. Returning origin instead\n.", x, y, z, t);
        return assign_position(0, 0, 0, 0);
    }

    return position;
}

void print_pos_vec(const pos_vec pos) {
    //	prints a position to the screen

    printf("x: %hu y: %hu z: %hu t: %hu\n", pos.i, pos.j, pos.k, pos.t);
}

inline pos_vec hop_pos_plus(const pos_vec u, const lorentz_idx mu) {
    //	Calculates the position immediately forward
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    #ifdef CHECK_POSITION_BOUNDS
        if(!position_mu_valid(u, mu)){
            printf("Position: ");
            print_pos_vec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }
    #endif

    pos_vec u_plus_muhat;

    unsigned short v;

    switch (mu) {
        case X_INDX:
            u_plus_muhat.i = ((v = u.i + 1) != N_SPC ? v : 0);
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = ((v = u.j + 1) != N_SPC ? v : 0);
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = ((v = u.k + 1) != N_SPC ? v : 0);
            u_plus_muhat.t = u.t;
            break;

        case T_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = ((v = u.t + 1) != N_T   ? v : 0);
            break;

    }

    return u_plus_muhat;
}

inline pos_vec hop_pos_minus(const pos_vec u, const lorentz_idx mu) {
    //	Calculates the position immediately behind
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    #ifdef CHECK_POSITION_BOUNDS
        if(!position_mu_valid(u, mu)){
            printf("Position: ");
            print_pos_vec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }
    #endif

    pos_vec u_minus_muhat;

    short v ;

     switch (mu) {
        case X_INDX:
            u_minus_muhat.i = ((v = u.i - 1) != -1 ? v : N_SPC - 1);
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = ((v = u.j - 1) != -1 ? v : N_SPC - 1);
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = ((v = u.k - 1) != -1 ? v : N_SPC - 1);   
            u_minus_muhat.t = u.t;
            break;

        case T_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = ((v = u.t - 1) != -1 ?   v : N_T - 1);
            break;
    }

    return u_minus_muhat;
}
