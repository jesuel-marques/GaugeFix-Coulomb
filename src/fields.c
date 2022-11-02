#include <stdio.h>
#include <stdlib.h>

#include <fields.h>
#include <flags.h>
#include <lattice.h>
#include <misc.h>
#include <SU3_ops.h>
#include <SU3_parameters.h>
#include <types.h>


extern short n_SPC;
extern short n_T;

extern int volume;


Mtrx3x3 *allocate_field(int elements, 
                        int size_of_elements){

    Mtrx3x3 *field = (Mtrx3x3 *) calloc(elements, size_of_elements);
	if(field == NULL){		
		return NULL;
	}

    return field;
}

void set_field_to_identity(Mtrx3x3 * restrict field, 
                           int elements){
    OMP_PARALLEL_FOR
    for(int i = 0; i < elements; i++)
        set_identity_3x3(field + i);
}

void copy_field(Mtrx3x3 * restrict field, 
                int elements, 
                Mtrx3x3 * restrict field_copy){
    OMP_PARALLEL_FOR
    for(int i = 0; i < elements; i++)
        copy_3x3(field + i, field_copy + i);
}

int reunitarize_field(Mtrx3x3 * restrict field,
                      int elements){
    int exit_status = 0;
    #pragma omp parallel for num_threads(NUM_THREADS) reduction (|:exit_status) schedule(dynamic)
    for(int i = 0; i < elements; i++)

        exit_status |= projection_SU3(field + i);

    if(exit_status != 0){
        return -1;
    }
    return 0;
}

double average_field_det(Mtrx3x3 *restrict field,
                         int elements) {
    double det = 0.0;
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) \
                        schedule(dynamic) reduction(+ : det)
        // Paralelizing by slicing the time extent
        for(int i = 0; i < elements; i++)

            det += determinant_3x3(field + i);

    det /= (DIM * volume);

    return det;
}


Mtrx3x3 *get_link(Mtrx3x3 *restrict U, 
                  const PosVec position, 
                  const LorentzIdx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    #ifdef CHECK_POSITION_BOUNDS
        if(position_mu_valid(position, mu))    
    #endif  //CHECK_POSITION_BOUNDS
            return U + GET_LINK_U(position, mu);
    #ifdef CHECK_POSITION_BOUNDS
        else{
            printf("Program reading config outside of allowed range.\n");
            exit(EXIT_FAILURE);
        }
    #endif  //CHECK_POSITION_BOUNDS
}

void get_link_matrix(      Mtrx3x3 *restrict U, 
                     const PosVec position, 
                     const LorentzIdx mu, 
                           Direction dir, 
                           Mtrx3x3 *restrict u) {
    // Gets forward or backward link at given position and mu
    // and copies it to u.

    #ifdef CHECK_POSITION_BOUNDS
    if(position_mu_valid(position, mu)){
    #endif  //CHECK_POSITION_BOUNDS
        if (dir == FRONT) {
        copy_3x3(get_link(U, position, mu), u);
        //	Link in the positive way is what is stored in U

        } else if (dir == REAR) {
        herm_conj_3x3(get_link(U, hop_pos_minus(position, mu), mu), u);
        //	U_(-mu)(n)=(U_mu(n-mu))^\dagger
        }
        else{
            printf("Direction is not valid.");
            exit(EXIT_FAILURE); 
        }
    #ifdef CHECK_POSITION_BOUNDS
    }
    else{
        printf("Program reading config outside of allowed range.\n");
        exit(EXIT_FAILURE);
    }
    #endif  //CHECK_POSITION_BOUNDS
    
}

Mtrx3x3 *get_gaugetransf(Mtrx3x3 *restrict G, 
                         const PosVec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position

    #ifdef CHECK_POSITION_BOUNDS
        if(position_valid(position)){
    #endif  //CHECK_POSITION_BOUNDS
            return G + GET_GT(position);
    #ifdef CHECK_POSITION_BOUNDS
        }
        else{
        
            printf("Program reading gauge-transformation outside of allowed position range.\n");
            exit(EXIT_FAILURE);
        }
    #endif  //CHECK_POSITION_BOUNDS
}