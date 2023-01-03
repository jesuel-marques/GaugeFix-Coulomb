#include <stdio.h>
#include <stdlib.h>

#include <fields.h>
#include <flags.h>
#include <lattice.h>
#include <misc.h>
#include <SU3_ops.h>

#include <types.h>


Mtrx3x3 *allocate_field(int elements, 
                        int size_of_elements){
    if(elements <= 0  || size_of_elements <= 0){
        return NULL;
    }
    
    const Mtrx3x3 *field = (Mtrx3x3 *) calloc(elements, size_of_elements);
	
    return field;
}


int set_field_to_identity(const Mtrx3x3 * restrict field, 
                           int elements){
    if(elements <= 0  || field == NULL){
        return -1;
    }

    OMP_PARALLEL_FOR
    for(int i = 0; i < elements; i++)
        set_identity_3x3(field + i);
    
    return 0;
}


int copy_field(const Mtrx3x3 * restrict field, 
                int elements, 
               const Mtrx3x3 * restrict field_copy){
    if(elements <= 0  || field == NULL || field_copy == NULL){
        return -1;
    }

    OMP_PARALLEL_FOR
    for(int i = 0; i < elements; i++)
        copy_3x3(field + i, field_copy + i);

    return 0;
}


int reunitarize_field(const Mtrx3x3 * restrict field,
                      int elements){
    if(elements <= 0  || field == NULL){
        return -1;
    }

    int exit_status = 0;
    #pragma omp parallel for num_threads(NUM_THREADS) \
            reduction (|:exit_status) schedule(dynamic)
    for(int i = 0; i < elements; i++)
        exit_status |= projection_SU3(field + i);

    if(exit_status != 0){
        return -1;
    }
    return 0;
}


double average_field_det(const Mtrx3x3 *restrict field,
                         int elements) {
    if(elements <= 0  || field == NULL){
        return -1;
    }
    double det = 0.0;

    #pragma omp parallel for num_threads(NUM_THREADS) \
            schedule(dynamic) reduction(+ : det)
        for(int i = 0; i < elements; i++)
            det += determinant_3x3(field + i);

    det /= (double)elements;

    return det;
}


Mtrx3x3 *get_link(const Mtrx3x3 *restrict U, 
                  const PosVec position, 
                  const LorentzIdx mu) {
    /* Does the pointer arithmetic to get a pointer to link at given position and mu */
    #ifdef CHECK_POSITION_BOUNDS
        if(position_mu_valid(position, mu))    
    #endif  //CHECK_POSITION_BOUNDS
            return GET_LINK(U, position, mu);
    #ifdef CHECK_POSITION_BOUNDS
        else{
            printf("Program reading config outside of allowed range.\n");
            exit(EXIT_FAILURE);
        }
    #endif  //CHECK_POSITION_BOUNDS
}


void get_link_matrix(const Mtrx3x3 *restrict U, 
                     const PosVec position, 
                     const LorentzIdx mu, 
                           Direction dir, 
                     const Mtrx3x3 *restrict u) {
    /* Gets forward or backward link at given position and mu
    and copies it to u. */
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


Mtrx3x3 *get_gaugetransf(const Mtrx3x3 *restrict G, 
                         const PosVec position) {
    /* Does the pointer arithmetic to get a pointer to 
       a gauge transformation at given position */
    #ifdef CHECK_POSITION_BOUNDS
        if(position_valid(position)){
    #endif  //CHECK_POSITION_BOUNDS
            return GET_GT(G, position);
    #ifdef CHECK_POSITION_BOUNDS
        }
        else{
        
            fprintf(stderr, "Program reading gauge-transformation"
                            "outside of allowed position range.\n");
            exit(EXIT_FAILURE);
            
        }
    #endif  //CHECK_POSITION_BOUNDS
}