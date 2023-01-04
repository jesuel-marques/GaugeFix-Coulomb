#include <stdio.h>
#include <stdbool.h>

#include <flags.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>

GeometricParameters lattice_param;

#define HOP_POS_SPA(v, u, i)    (v).i = ((u).i != lattice_param.n_SPC - 1 ? \
                                                   (u).i + 1 : 0);
#define HOP_POS_TIME(v, u)      (v).t  = ((u).t  != lattice_param.n_T   - 1 ? \
                                                   (u).t  + 1 : 0);

#define HOP_MIN_SPA(v, u, i)    (v).i = ((u).i != 0 ? \
                                         (u).i - 1 : lattice_param.n_SPC - 1);
#define HOP_MIN_TIME(v, u)      (v).t  = ((u).t  != 0 ? \
                                         (u).t  - 1 : lattice_param.n_T - 1);


int initGeometricParameters(const short n_s, const short n_t) {

    lattice_param.n_SPC 		= n_s;
	lattice_param.n_T   		= n_t;

	if(lattice_param.n_SPC <= 0 || lattice_param.n_T <= 0) {
		lattice_param.error = 1;
	}

	lattice_param.spatial_volume = lattice_param.n_SPC * 
                                   lattice_param.n_SPC * 
                                   lattice_param.n_SPC;	
	lattice_param.volume 	     = lattice_param.spatial_volume * 
                                   lattice_param.n_T;

	lattice_param.amount_of_links  = DIM * lattice_param.volume;
	lattice_param.amount_of_points = lattice_param.volume;

    if(validGeometricParametersQ()){
        return 0;
    }
    else{
        return 1;
    }
}

bool validGeometricParametersQ(void){
    short ns = lattice_param.n_SPC;
    short nt = lattice_param.n_T;

    if(lattice_param.n_SPC <= 0                                 || 
       lattice_param.n_T   <= 0                                 ||
       lattice_param.error == 1                                 ||
       lattice_param.spatial_volume != ns * ns * ns             ||
       lattice_param.volume != ns * ns * ns * nt                ||
       lattice_param.amount_of_links != DIM * ns * ns * ns * nt ||
       lattice_param.amount_of_points != ns * ns * ns * nt        ) {
        
        lattice_param.error = 1;
        return false;

    }
    else{

        return true;

    }
	
}

inline bool validPositionQ(PosVec position) {
    if (position.t < lattice_param.n_T   &&
        position.t >= 0                  &&
        position.i < lattice_param.n_SPC &&
        position.i >= 0                  &&
        position.j < lattice_param.n_SPC &&
        position.j >= 0                  &&
        position.k < lattice_param.n_SPC &&
        position.k >= 0                    ) {
        return true;
    }
    else{

        return false;

    }
}


inline bool positionmuValidQ(PosVec position, 
                              LorentzIdx mu) {
    if (validPositionQ(position) &&
        (mu == T_INDX || mu == X_INDX || 
         mu == Y_INDX || mu == Z_INDX   )) {

        return true;

    }

    return false;
}


int mod(short a, unsigned short b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}


PosVec makePositionValid(PosVec position) {
    PosVec valid_position;
    
    if(!validGeometricParametersQ()){
        fprintf(stderr, "Not possible to make position valid." 
                        "Error in lattice parameters\n"
                        "Returning same possibly invalid position.");
        return position;
    }

    valid_position.i = mod(position.i, lattice_param.n_SPC);
    valid_position.j = mod(position.j, lattice_param.n_SPC);
    valid_position.k = mod(position.k, lattice_param.n_SPC);
    valid_position.t = mod(position.t, lattice_param.n_T  );

    return valid_position;   
}


PosVec assignPosition(const PosIndex x, 
                      const PosIndex z, 
                      const PosIndex y, 
                      const PosIndex t) {
    /* assigns x, y, z and t to a position vector */
    PosVec position = {.i = x, .j = y, .k = z, .t = t};
    
    if(!validPositionQ(position)) {
        position = makePositionValid(position);
    }   /* if position not valid, make it valid and assign the valid one */
    return position;
}

void printPosVec(const PosVec pos) {
    /* prints a position to the screen */
    printf("x: %hu y: %hu z: %hu t: %hu\n", pos.i, pos.j, pos.k, pos.t);
}


inline PosVec hopPosPlus(const PosVec u, 
                         const LorentzIdx mu) {
    /* 	Calculates the position immediately forward
    	in the direction mu, taken into account the
    	periodic boundary conditions. */
    #ifdef CHECK_POSITION_BOUNDS
        if(!positionmuValidQ(u, mu)) {
            printf("Position: ");
            printPosVec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }

        if(!validGeometricParametersQ()){
            fprintf(stderr, "Error in lattice parameters\n");
            exit(EXIT_FAILURE);
        }
    #endif  //CHECK_POSITION_BOUNDS

    PosVec u_plus_muhat = u;

    switch (mu) {
        case X_INDX:
            HOP_POS_SPA(u_plus_muhat, u, i);
            break;

        case Y_INDX:
            HOP_POS_SPA(u_plus_muhat, u, j);
            break;

        case Z_INDX:
            HOP_POS_SPA(u_plus_muhat, u, k);
            break;

        case T_INDX:
            HOP_POS_TIME(u_plus_muhat, u);
            break;

    }

    return u_plus_muhat;
}


inline PosVec hopPosMinus(const PosVec u, 
                          const LorentzIdx mu) {
    /* Calculates the position immediately behind
    in the direction mu, taken into account the
    periodic boundary conditions. */
    #ifdef CHECK_POSITION_BOUNDS
        if(!positionmuValidQ(u, mu)) {
            printf("Position: ");
            printPosVec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);

            if(!validGeometricParametersQ()){
                fprintf(stderr, "Error in lattice parameters\n");
                exit(EXIT_FAILURE);
            }
        }
    #endif  //CHECK_POSITION_BOUNDS

    PosVec u_minus_muhat = u;

    switch (mu) {
        case X_INDX:
            HOP_MIN_SPA(u_minus_muhat, u, i);
            break;

        case Y_INDX:
            HOP_MIN_SPA(u_minus_muhat, u, j);
            break;

        case Z_INDX:
            HOP_MIN_SPA(u_minus_muhat, u, k);
            break;

        case T_INDX:
            HOP_MIN_TIME(u_minus_muhat, u);
            break;

    }

    return u_minus_muhat;
}