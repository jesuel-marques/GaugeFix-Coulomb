/*
    fields_io takes care of input and output of lattice fields.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br
    
 */

#include <SU2_ops.h>
#include <SU3_ops.h>
#include <convert_precision_io.h>
#include <fields.h>
#include <fields_io.h>
#include <../settings.h>
#include <geometry.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

extern GeometricParameters lattice_param;

/*============================JONIVAR'S CODE===============================*/
/* conversion between little and big endian */

void swap_block(int *buffer, size_t length) {
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

void swap_block_double(double *buffer, size_t length) {
    union swapper {
        double double_number;
        char pos[8];
    } a, b;

    for (size_t i = 0; i < length; i++) {
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

int byteSwap(void *strip, size_t size, size_t length) {
    switch (size) {
        case sizeof(float): /* = 4 */
            swap_block(strip, length / size);
            break;
        case sizeof(double): /* = 8 */
            swap_block_double(strip, length / size);
            break;
        case 1:
            break;
        default: /* ERROR: unknown size!! */
            return -1;
    }

    return 0;
}
/*=============================================================================*/

/* Gets a pointer to a link in input precision at the given position and mu. */
InCfgMtrx *getLinkIn(InCfgMtrx *U_in, const PosVec position, const LorentzIdx mu) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * GET_LINK.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U_in:      SU(3) gluon field in input precision,
     * PosVec position:     position at which the link variable is requested,
     * LorentzIdx mu:       Lorentz direction for which the link variable is requested.
     *
     * Returns:
     * =======
     * The gauge link at that particular point and mu as a pointer to a 3x3 matrix in
     * input precision.
     *
     */

    return GET_LINK(U_in, position, mu);
}

/* Gets a pointer to a link in output precision at the given position and mu. */
OutCfgMtrx *getLinkOut(OutCfgMtrx *U_out, const PosVec position, const LorentzIdx mu) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * GET_LINK.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * U_out:     SU(3) gluon field in output precision,
     * PosVec position:     position at which the link variable is requested,
     * LorentzIdx mu:       Lorentz direction for which the link variable is requested.
     *
     * Returns:
     * =======
     * The gauge link at that particular point and mu as a pointer to a 3x3 matrix in
     * output precision.
     *
     */

    return GET_LINK(U_out, position, mu);
}

/* Gets a pointer to the gauge transformation in input precision at the given position. 
*/
InGTMtrx *getGaugetransfIn(InGTMtrx *restrict G_in, const PosVec position) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * GET_GT.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * G_in:      SU(3) gauge transformation field in input precision,
     * PosVec position:     position at which the gauge transformation is requested.
     *
     * Returns:
     * =======
     * The gauge transformation at that particular point as a pointer to a 3x3 matrix in
     * input precision.
     *
     */

    return GET_GT(G_in, position);
}

/* Gets a pointer to the gauge transformation in output precision at the given
 position. */
OutGTMtrx *getGaugetransfOut(OutGTMtrx *restrict G_out, const PosVec position) {
    /*
     * Calls:
     * =====
     *
     * Macros:
     * ======
     * GET_GT.
     *
     * Global Variables:
     * ================
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * G_out:     SU(3) gauge transformation field in output precision,
     * PosVec position:     position at which the gauge transformation is requested.
     *
     * Returns:
     * =======
     * The gauge transformation at that particular point as a pointer to a 3x3 matrix in
     * output precision.
     *
     */

    return GET_GT(G_out, position);
}

/* Reads a file with filename config_filename and loads a gluon-field from it. */
int loadConfig(const Mtrx3x3 *restrict U, char * config_filename) {

    /*
     * Calls:
     * =====
     * fopen, fprintf, fread, fclose, feof, fgetc,
     * if CONV_CFG_TO_WORKING_PRECISION set
     * calloc, free,
     * byteSwap, 
     * convertCfg_in_work.
     *
     * Macros:
     * ======
     * CONV_CFG_TO_WORKING_PRECISION, TEST_ALLOCATION, NEED_BYTE_SWAP_IN, DIM.
     *
     * Global Variables:
     * ================
     * lattice_param
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * restrict U:    the gluon field to be written,
     * char * config_filename:  the filename of the file which will store the
     *                          gluon field.
     *
     * Returns:
     * =======
     * 0 on success, negative integer on failure.
     */
    
    FILE *config_file;

    if((config_file = fopen(config_filename, "rb")) == NULL) {
        fprintf(stderr, "Error: Problem opening config file %s.\n", config_filename);
        return -1;
    }

    const InCfgMtrx *U_in;

    #ifdef CONV_CFG_TO_WORKING_PRECISION
        U_in = (InCfgMtrx *)calloc(lattice_param.volume * DIM, sizeof(InCfgMtrx));
        if(U_in == NULL) {
            return -3;
        }
    #else   // CONV_CFG_TO_WORKING_PRECISION
        U_in = (InCfgMtrx *)U;
    #endif  // CONV_CFG_TO_WORKING_PRECISION

        if(fread(U_in, sizeof(InCfgMtrx), lattice_param.volume * DIM, config_file) !=
                                          lattice_param.volume * DIM) {
            fprintf(stderr, "Error: Reading from file %s.\n", config_filename);

    #ifdef CONV_CFG_TO_WORKING_PRECISION
            free(U_in);
    #endif  // CONV_CFG_TO_WORKING_PRECISION
        fclose(config_file);
        return -2;
    }
    fgetc(config_file);

    if(!feof(config_file)) {
        fprintf(stderr,
                "Error: File has not been read till the end.\n Check lattice sizes.\n");
    #ifdef CONV_CFG_TO_WORKING_PRECISION
        free(U_in);
    #endif  // CONV_CFG_TO_WORKING_PRECISION
        fclose(config_file);
        return -3;
    }

    fclose(config_file);

    #ifdef NEED_BYTE_SWAP_IN
        if(byteSwap(U_in, sizeof(InScalar) / 2, 
                    lattice_param.volume * DIM * sizeof(InCfgMtrx))) {
            fprintf(stderr, "Error: Problem with the byteSwap. \n");
    #ifdef CONV_CFG_TO_WORKING_PRECISION
            free(U_in);
    #endif  // CONV_CFG_TO_WORKING_PRECISION
            return -4;
        }
    #endif  // NEED_BYTE_SWAP_IN

    #ifdef CONV_CFG_TO_WORKING_PRECISION
        convertCfg_in_work(U_in, U);
        free(U_in);
    #endif  // CONV_CFG_TO_WORKING_PRECISION

    return 0;
}

/* Reads a file with filename gauge_transf_filename and loads a gauge transformation
   from it. */
int loadGaugeTransf(Mtrx3x3 *restrict G, char *gauge_transf_filename) {

    /*
     * Calls:
     * =====
     * fopen, fprintf, fread, fclose, feof, fgetc,
     * if CONV_GT_TO_WORKING_PRECISION set
     * calloc, free,
     * convertGT_in_work.
     *
     * Macros:
     * ======
     * CONV_GT_TO_WORKING_PRECISION,
     *
     * if CONV_GT_TO_WORKING_PRECISION set
     * TEST_ALLOCATION.
     *
     * Global Variables:
     * ================
     * lattice_param
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * restrict G:            the gauge transformation field to be read into,
     * char * gauge_transf_filename:    the filename of the file which stores the
     *                                  gauge transformation.
     *
     * Returns:
     * =======
     * 0 on success, negative integer on failure.
     */
    
    const InGTMtrx *G_in;

    #ifdef CONV_GT_TO_WORKING_PRECISION
        G_in = (InGTMtrx *)calloc(lattice_param.volume, sizeof(InGTMtrx));
        if(TEST_ALLOCATION(G_in)) {
            return -3;
        }
    #else   // CONV_GT_TO_WORKING_PRECISION
        G_in = (InGTMtrx *)G;
    #endif  // CONV_GT_TO_WORKING_PRECISION

    const FILE *gaugetransf_file;

    if((gaugetransf_file = fopen(gauge_transf_filename, "rb")) == NULL) {
        fprintf(stderr, 
                "Error: Problem opening file %s to load gauge transformation.\n",
                gauge_transf_filename);

    #ifdef CONV_GT_TO_WORKING_PRECISION
        free(G_in);
    #endif  // CONV_GT_TO_WORKING_PRECISION

        return -1;
    }

    if(fread(G_in, sizeof(InGTMtrx), lattice_param.volume, gaugetransf_file) !=
                                     lattice_param.volume) {
        fprintf(stderr, 
                "Error: Problem reading file %s to load gauge transformation.\n",
                gauge_transf_filename);

    #ifdef CONV_GT_TO_WORKING_PRECISION
        free(G_in);
    #endif  // CONV_GT_TO_WORKING_PRECISION

        fclose(gaugetransf_file);
        return -2;
    }

    fgetc(gaugetransf_file);

    if(!feof(gaugetransf_file)) {
        fprintf(stderr,
                "Error: File has not been read till the end.\nCheck lattice sizes.\n");

    #ifdef CONV_GT_TO_WORKING_PRECISION
        free(G_in);
    #endif  // CONV_GT_TO_WORKING_PRECISION

        fclose(gaugetransf_file);
        return -2;
    }

    fclose(gaugetransf_file);

    #ifdef CONV_GT_TO_WORKING_PRECISION
        convertGT_in_work(G_in, G);
        free(G_in);
    #endif  // CONV_GT_TO_WORKING_PRECISION

    return 0;
}

/* Creates a file with filename config_filename and stores a gluon-field in it. */
int writeConfig(Mtrx3x3 *restrict U, char *config_filename) {
    /*
     * Calls:
     * =====
     * fopen, fprintf, fwrite, fclose,
     * if CONV_CFG_FROM_WORKING_PRECISION set
     * calloc, free,
     * byteSwap, 
     * convertCfg_work_out.
     *
     * Macros:
     * ======
     * CONV_CFG_FROM_WORKING_PRECISION, TEST_ALLOCATION, NEED_BYTE_SWAP_OUT, DIM.
     *
     * Global Variables:
     * ================
     * lattice_param
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * restrict U:    the gluon field to be written,
     * char * config_filename:  the filename of the file which will store the
     *                          gluon field.
     *
     * Returns:
     * =======
     * 0 on success, negative integer on failure.
     */

    const OutCfgMtrx *U_out;

#ifdef CONV_CFG_FROM_WORKING_PRECISION

    U_out = (OutCfgMtrx *)calloc(lattice_param.volume * DIM, sizeof(OutCfgMtrx));
    TEST_ALLOCATION(U_out);

    convertCfg_work_out(U, U_out);

#else   // CONV_CFG_FROM_WORKING_PRECISION
    U_out = (OutCfgMtrx *)U;
#endif  // CONV_CFG_FROM_WORKING_PRECISION

#ifdef NEED_BYTE_SWAP_OUT
    if(byteSwap(U_out, sizeof(OutScalar) / 2,
                 lattice_param.volume * DIM * sizeof(OutCfgMtrx))) {
        fprintf(stderr, "Error: Problem with the byte swap. Unknown size.\n");
#ifdef CONV_CFG_FROM_WORKING_PRECISION
        free(U_out);
#endif  // CONV_CFG_FROM_WORKING_PRECISION
        return -3;
    }
#endif  // NEED_BYTE_SWAP_OUT

    const FILE *config_file;

    if((config_file = fopen(config_filename, "wb")) == NULL) {
        fprintf(stderr, "Error: Problem creating file %s for config.\n", 
                        config_filename);

#ifdef CONV_CFG_FROM_WORKING_PRECISION
        free(U_out);
#endif  // CONV_CFG_FROM_WORKING_PRECISION

        return -1;
    }

    if(fwrite(U_out, sizeof(OutCfgMtrx), lattice_param.volume * DIM, config_file) !=
        lattice_param.volume * DIM) {
        fclose(config_file);

#ifdef CONV_CFG_FROM_WORKING_PRECISION
        free(U_out);
#endif  // CONV_CFG_FROM_WORKING_PRECISION

        return -2;
    }

    fclose(config_file);

#ifdef CONV_CFG_FROM_WORKING_PRECISION
    free(U_out);
#endif  // CONV_CFG_FROM_WORKING_PRECISION

    return 0;
}

/* Creates a file with filename gauge_transf_filename and stores a gauge transformation
   in it. */
int writeGaugeTransf(Mtrx3x3 *restrict G, char *gauge_transf_filename) {

    /*
     * Calls:
     * =====
     * fopen, fprintf, fwrite, fclose,
     * if CONV_GT_FROM_WORKING_PRECISION set
     * calloc, free,
     * convertGT_work_out.
     *
     * Macros:
     * ======
     * CONV_GT_FROM_WORKING_PRECISION,
     *
     * if CONV_GT_FROM_WORKING_PRECISION set
     *  TEST_ALLOCATION.
     *
     * Global Variables:
     * ================
     * lattice_param
     *
     * Parameters:
     * ==========
     * Mtrx3x3 * restrict G:            the gauge transformation to be written,
     * char * gauge_transf_filename:    the filename of the file which will store the
     *                                  gauge transformation.
     *
     * Returns:
     * =======
     * 0 on success, negative integer on failure.
     */

    const OutGTMtrx *G_out;

    #ifdef CONV_GT_FROM_WORKING_PRECISION
        G_out = (OutGTMtrx *)calloc(lattice_param.volume, sizeof(OutGTMtrx));
        if(TEST_ALLOCATION(G_out)) {
            return -3;
        }
        convertGT_work_out(G, G_out);
    #else   // CONV_GT_FROM_WORKING_PRECISION
        G_out = (OutGTMtrx *)G;
    #endif  // CONV_GT_FROM_WORKING_PRECISION

    const FILE *gaugetransf_file;

    if((gaugetransf_file = fopen(gauge_transf_filename, "wb")) == NULL) {
        #ifdef CONV_GT_FROM_WORKING_PRECISION
            free(G_out);
        #endif  // CONV_GT_FROM_WORKING_PRECISION

        return -1;
    }

    if(fwrite(G_out, sizeof(OutGTMtrx), lattice_param.volume, gaugetransf_file) !=
        lattice_param.volume) {
        fclose(gaugetransf_file);

        #ifdef CONV_GT_FROM_WORKING_PRECISION
            free(G_out);
        #endif  // CONV_GT_FROM_WORKING_PRECISION

        return -2;
    }

    fclose(gaugetransf_file);

    #ifdef CONV_GT_FROM_WORKING_PRECISION
        free(G_out);
    #endif  // CONV_GT_FROM_WORKING_PRECISION

    return 0;
}