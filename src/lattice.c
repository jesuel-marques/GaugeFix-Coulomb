#include <stdio.h>
#include <stdbool.h>

#include <flags.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>

GeometricParameters lattice_param;

PosVec *neighbour_table;

static void makeNeighbourTable(){

    neighbour_table = (PosVec *) malloc(lattice_param.amount_of_points
                                        * DIM * 2 * sizeof(PosVec));
    short t;
    LOOP_TEMPORAL_PARALLEL(t){
        PosVec position = {.pos={0, 0, 0, t}};
        PosVec neighbour = {.pos={0, 0, 0, 0}};

        LorentzIdx mu;
        LOOP_SPATIAL(position){
            LOOP_LORENTZ(mu){
                neighbour = position;
                neighbour.pos[mu]--;

                *(neighbour_table + ((((position.pos[T_INDX]  * lattice_param.n_SPC 
                                      + position.pos[Z_INDX]) * lattice_param.n_SPC
                                      + position.pos[Y_INDX]) * lattice_param.n_SPC
                                      + position.pos[X_INDX]) * DIM 
                                      + mu) * 2 + REAR) = makePeriodicBound(neighbour);
            
                neighbour = position;
                neighbour.pos[mu]++;

                *(neighbour_table + ((((position.pos[T_INDX]  * lattice_param.n_SPC 
                                      + position.pos[Z_INDX]) * lattice_param.n_SPC
                                      + position.pos[Y_INDX]) * lattice_param.n_SPC
                                      + position.pos[X_INDX]) * DIM 
                                      + mu) * 2 + FRONT) = makePeriodicBound(neighbour);
            }
        }
    }

}

void finalizeGeometry(){
    free(neighbour_table);
}

int initGeometry(const short n_s, const short n_t) {

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

    if(!validGeometricParametersQ()){
        return 1;
    }

    makeNeighbourTable();

    return 0;
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

    for(LorentzIdx mu = X_INDX; mu < DIM; mu++){
        if(position.pos[mu] >= (mu != T_INDX ? 
                               lattice_param.n_SPC : 
                               lattice_param.n_T) &&
           position.pos[mu] < 0 )
            return false;
    }

    return true;
}


inline bool positionmuValidQ(PosVec position, 
                             LorentzIdx mu) {
    if(validPositionQ(position) &&
        (mu == T_INDX || mu == X_INDX || 
         mu == Y_INDX || mu == Z_INDX   )) {

        return true;

    }

    return false;
}


inline int mod(short a, unsigned short b)
{
    short r = a % b;
    return r < 0 ? r + b : r;
}


PosVec makePeriodicBound(PosVec position) {
    PosVec valid_position;
    
    // if(!validGeometricParametersQ()){
    //     fprintf(stderr, "Error in geometric parameters\n";
    //     return position;
    // }

    for(LorentzIdx mu = X_INDX; mu < DIM; mu++){
        valid_position.pos[mu] = mod(position.pos[mu], mu != T_INDX ? 
                                                       lattice_param.n_SPC :
                                                       lattice_param.n_T);
    }
    
    return valid_position;   
}


PosVec assignPosition(const PosIndex x, 
                      const PosIndex z, 
                      const PosIndex y, 
                      const PosIndex t) {
    /* assigns x, y, z and t to a position vector */
    PosVec position = {.pos = {x, y, z, t}};
    
    if(!validPositionQ(position)) {
        position = makePeriodicBound(position);
    }   /* if position not valid, make it valid and assign the valid one */
    return position;
}

void printPosVec(const PosVec position) {
    /* prints a position to the screen */
    printf("x: %hu y: %hu z: %hu t: %hu\n", position.pos[X_INDX],
                                            position.pos[Y_INDX],
                                            position.pos[Z_INDX],
                                            position.pos[T_INDX] );

}


inline PosVec getNeighbour(PosVec position, LorentzIdx mu, Direction dir){
    return *(neighbour_table + ((((position.pos[T_INDX]  * lattice_param.n_SPC 
                                  +position.pos[Z_INDX]) * lattice_param.n_SPC
                                  +position.pos[Y_INDX]) * lattice_param.n_SPC
                                  +position.pos[X_INDX]) * DIM + mu) * 2 + dir);
}