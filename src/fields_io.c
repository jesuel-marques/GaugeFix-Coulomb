#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <convert_precision_io.h>
#include <fields.h>
#include <fields_io.h>
#include <flags.h>
#include <lattice.h>
#include <misc.h>
#include <SU2_ops.h>
#include <SU3_ops.h>
#include <SU3_parameters.h>


extern short n_SPC;
extern short n_T;

extern int volume;

/*============================JONIVAR'S CODE===============================*/

void block_swap(int *buffer, size_t length) {
    size_t i;
    union swapper {
        int integer;
        char pos[4];
    } a, b;

    for (i = 0; i < length; i++) {
        a.integer = *buffer;
        b.pos[0] = a.pos[3];
        b.pos[1] = a.pos[2];
        b.pos[2] = a.pos[1];
        b.pos[3] = a.pos[0];
        *buffer = b.integer;
        buffer++;
    }
}


void block_swap_double(double *buffer, size_t length) {
    size_t i;
    union swapper {
        double double_number;
        char pos[8];
    } a, b;

    for (i = 0; i < length; i++) {
        a.double_number = *buffer;
        b.pos[0] = a.pos[7];
        b.pos[1] = a.pos[6];
        b.pos[2] = a.pos[5];
        b.pos[3] = a.pos[4];
        b.pos[4] = a.pos[3];
        b.pos[5] = a.pos[2];
        b.pos[6] = a.pos[1];
        b.pos[7] = a.pos[0];
        *buffer = b.double_number;
        buffer++;
    }
}


int byte_swap(void *strip, size_t size, size_t length) {
    switch (size) {
        case sizeof(float): /* = 4 */
            block_swap(strip, length / size);
            break;
        case sizeof(double): /* = 8 */
            block_swap_double(strip, length / size);
            break;
        case 1:
            break;
        default: /* ERROR: unknown size!! */
            return -1;
    }

    return 0;
}
/*==========================================================================*/


InCfgMtrx *get_link_in(InCfgMtrx *U, 
                       const PosVec position,
                       const LorentzIdx mu) {
    //	Does the pointer arithmetic to get a pointer
    //  to link at given position and mu
    return U + GET_LINK_U(position, mu);
}


OutCfgMtrx *get_link_out(OutCfgMtrx *U, 
                         const PosVec position,
                         const LorentzIdx mu) {
    //	Does the pointer arithmetic to get a pointer
    //  to link at given position and mu
    return U + GET_LINK_U(position, mu);
}


InGTMtrx  *get_gaugetransf_in (InGTMtrx  * restrict G_in,
                               const PosVec position) {
    //	Does the pointer arithmetic to get a pointer 
    //  to a gauge-transformation at given position
    return G_in + GET_GT(position);
}


OutGTMtrx *get_gaugetransf_out(OutGTMtrx * restrict G_out,
                               const PosVec position) {
    //	Does the pointer arithmetic to get a pointer 
    //  to a gauge-transformation at given position
    return G_out + GET_GT(position);
}


int SU3_load_config(Mtrx3x3 * restrict U, 
                    char * config_filename) {
    //	Loads a link configuration from the file with filename to U.
    FILE *config_file;

    if ((config_file = fopen(config_filename, "rb")) == NULL) {

        fprintf(stderr, "Error: Problem opening config file %s.\n", config_filename);

        return -1;
    }

    InCfgMtrx *U_in;

    #ifdef CONV_CFG_TO_WORKING_PRECISION

        U_in = (InCfgMtrx *)calloc(volume * DIM, sizeof(InCfgMtrx));
        if(TEST_ALLOCATION(U_in)){

            return -1;
        }

    #else   //CONV_CFG_TO_WORKING_PRECISION
        U_in = (InCfgMtrx *)U;
    #endif  //CONV_CFG_TO_WORKING_PRECISION

    if (fread(U_in, sizeof(InCfgMtrx), volume * DIM, config_file) != volume * DIM) {

        fprintf(stderr, "Error: Reading from file %s.\n", config_filename);
        #ifdef CONV_CFG_TO_WORKING_PRECISION 
             free(U_in);
        #endif  //CONV_CFG_TO_WORKING_PRECISION
        fclose(config_file);
        return -1;

    }
    fgetc(config_file);

    if (!feof(config_file)) {
        fprintf(stderr, "Error: File has not been read till the end.\n"
                        "Check lattice sizes.\n");
        #ifdef CONV_CFG_TO_WORKING_PRECISION 
             free(U_in);
        #endif  //CONV_CFG_TO_WORKING_PRECISION
        fclose(config_file);
        return -1;
    }

    fclose(config_file);

    #ifdef NEED_BYTE_SWAP_IN
        if(byte_swap(U_in, sizeof(InScalar) / 2, volume * DIM * sizeof(InCfgMtrx))){
            fprintf(stderr, "Error: Problem with the byte_swap. Unknown size.\n");
            #ifdef CONV_CFG_TO_WORKING_PRECISION
                free(U_in);
            #endif  //CONV_CFG_TO_WORKING_PRECISION
            return -1;
        }
    #endif  //NEED_BYTE_SWAP_IN

    #ifdef CONV_CFG_TO_WORKING_PRECISION
        SU3_convert_cfg_in_work(U_in, U);
        free(U_in);
    #endif  //CONV_CFG_TO_WORKING_PRECISION

    return 0;
}


int SU3_load_gauge_transf(Mtrx3x3 * restrict G, 
                          char    * gauge_transf_filename) {
    //	Loads a gauge transformation to G.
    InGTMtrx *G_in;

    #ifdef CONV_GT_TO_WORKING_PRECISION
        G_in = (InGTMtrx *)calloc( volume , sizeof(InGTMtrx));
        if(TEST_ALLOCATION(G_in)){
            return -1;
        }
    #else   //CONV_GT_TO_WORKING_PRECISION
        G_in = (InGTMtrx *)G;
    #endif  //CONV_GT_TO_WORKING_PRECISION
    
    FILE *gaugetransf_file;

    printf("Loading: %s.\n", gauge_transf_filename);
    if ((gaugetransf_file = fopen(gauge_transf_filename, "rb")) == NULL) {

        fprintf(stderr, "Error: Problem opening file %s"
                        "to load gauge transformation.\n", gauge_transf_filename);

        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif  //CONV_GT_TO_WORKING_PRECISION

        return -1;

    }

    if (fread(G_in, sizeof(InGTMtrx), volume, gaugetransf_file) != volume ) {

        fprintf(stderr, "Error: Problem reading file %s"
                        "to load gauge transformation.\n", gauge_transf_filename);

        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif  //CONV_GT_TO_WORKING_PRECISION

        fclose(gaugetransf_file);
        return -1;

    }

    fgetc(gaugetransf_file);

    if (!feof(gaugetransf_file)) {

        fprintf(stderr, "Error: File has not been read till the end."
                        "Check lattice sizes.\n");
        
        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif  //CONV_GT_TO_WORKING_PRECISION

        fclose(gaugetransf_file);
        return -1;
    }

    fclose(gaugetransf_file);

    #ifdef CONV_GT_TO_WORKING_PRECISION
        SU3_convert_gt_in_work(G_in, G);
        free(G_in);
    #endif  //CONV_GT_TO_WORKING_PRECISION

    return 0;
}


int SU3_write_config(Mtrx3x3 * restrict U, 
                     char * config_filename) {
    //  Loads a link configuration from the file with filename to U.
    OutCfgMtrx *U_out;

    #ifdef CONV_CFG_FROM_WORKING_PRECISION

        U_out = (OutCfgMtrx *)calloc(volume * DIM, sizeof(OutCfgMtrx));
        TEST_ALLOCATION(U_out);

        SU3_convert_cfg_work_out(U, U_out);

    #else   //CONV_CFG_FROM_WORKING_PRECISION
        U_out = (OutCfgMtrx *)U;
    #endif  //CONV_CFG_FROM_WORKING_PRECISION

    #ifdef NEED_BYTE_SWAP_OUT
        if (byte_swap(U_out, sizeof(OutScalar) / 2, volume * DIM * sizeof(OutCfgMtrx))){
            fprintf(stderr, "Error: Problem with the byte_swap. Unknown size.\n");
            #ifdef CONV_CFG_FROM_WORKING_PRECISION
                free(U_out);
            #endif  //CONV_CFG_FROM_WORKING_PRECISION
            return -1;
        }
    #endif  //NEED_BYTE_SWAP_OUT

    FILE *config_file;

    printf("Creating: %s.\n", config_filename);
    if ((config_file = fopen(config_filename, "wb")) == NULL) {

        fprintf(stderr, "Error: Problem creating file %s"
                        "for config.\n", config_filename);

        #ifdef CONV_CFG_FROM_WORKING_PRECISION
            free(U_out);
        #endif  //CONV_CFG_FROM_WORKING_PRECISION

        return -1;

    }

    if (fwrite(U_out, sizeof(OutCfgMtrx), volume * DIM, config_file) != volume * DIM) {

        fclose(config_file);

        #ifdef CONV_CFG_FROM_WORKING_PRECISION
            free(U_out);
        #endif  //CONV_CFG_FROM_WORKING_PRECISION

        return -1;
    }

    fclose(config_file);

    #ifdef CONV_CFG_FROM_WORKING_PRECISION
        free(U_out);
    #endif  //CONV_CFG_FROM_WORKING_PRECISION

    return 0;
}


int SU3_write_gauge_transf(Mtrx3x3 * restrict G, 
                           char * gauge_transf_filename) {
    //  Loads a link configuration from the file with filename to U.
    OutGTMtrx *G_out;

    #ifdef CONV_GT_FROM_WORKING_PRECISION
        G_out = (OutGTMtrx *)calloc(volume, sizeof(OutGTMtrx));
        if(TEST_ALLOCATION(G_out)){
            
            return -1;
        
        }
        SU3_convert_gt_work_out(G, G_out);
    #else   //CONV_GT_FROM_WORKING_PRECISION
        G_out = (OutGTMtrx *)G;
    #endif  //CONV_GT_FROM_WORKING_PRECISION

    FILE *gaugetransf_file;

    if ((gaugetransf_file = fopen(gauge_transf_filename, "wb")) == NULL) {

        fprintf(stderr, "Error: Problem opening file %s \
                         to store gauge transformation.\n", gauge_transf_filename);
        
        #ifdef CONV_GT_FROM_WORKING_PRECISION
            free(G_out);
        #endif  //CONV_GT_FROM_WORKING_PRECISION

        return -1;
    }

    if (fwrite(G_out, sizeof(OutGTMtrx), volume, gaugetransf_file) != volume ) {

        fprintf(stderr, "Error: Problem writing gauge transformation to file %s.\n",
                gauge_transf_filename);

        fclose(gaugetransf_file);

        #ifdef CONV_GT_FROM_WORKING_PRECISION

            free(G_out);
            
        #endif  //CONV_GT_FROM_WORKING_PRECISION

        return -1;
    }

    fclose(gaugetransf_file);

    #ifdef CONV_GT_FROM_WORKING_PRECISION

        free(G_out);
        
    #endif  //CONV_GT_FROM_WORKING_PRECISION

    return 0;
}