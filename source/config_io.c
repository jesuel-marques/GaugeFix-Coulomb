
#include <complex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <string.h>



#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters
#include "../SU3_parameters.h"              //	Simulation parameters
#include "SU2_ops.h"                        //	SU(2) operations
#include "SU3_ops.h"                        //	SU(3) operations
                                            //	positions and links on the lattice.

#include "lattice.h"  //	Initialization functions and calculations of

extern char config_template[];
extern char configs_dir_name_in[];
extern char configs_dir_name_out[];
extern char extension_config_in[];
extern char extension_config_out[];

extern char extension_gt_in[];
extern char extension_gt_out[];

extern const int config_exception_list[];



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

void greeter_function(char *program_name) {
    printf("Program %s compiled at %s on %s\n", program_name, __TIME__, __DATE__);
    printf("Using C version: %ld\n", __STDC_VERSION__);
}

void handle_input(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Input names for input and output config directory and template");

        if (argc > 4) {
            fprintf(stderr, " only");
        }
    }

    if (argc != 4) {
        fprintf(stderr, ".\n");
        exit(EXIT_FAILURE);
    }

    char size_directory[10];

    sprintf(size_directory, "%dx%d/", N_SPC, N_T);

    strcpy(configs_dir_name_in, argv[1]);
    strcpy(configs_dir_name_out, argv[2]);

    strcat(configs_dir_name_in, size_directory);
    strcat(configs_dir_name_out, size_directory);

    strcpy(config_template, argv[3]);
}

bool is_in_exception_list(const int config_nr) {
    for (int i = 0; config_exception_list[i] != -1; i++) {
        if (config_nr == config_exception_list[i]) {
            return true;
        }
    }

    return false;
}

void create_output_directory(void) {
    char command_create_directory[MAX_LENGTH_NAME];
    sprintf(command_create_directory, "mkdir %s", configs_dir_name_out);

    if (system(command_create_directory)) {
        fprintf(stderr, "\nError creating directory %s to store gauge fixed configs\n", configs_dir_name_out);
        exit(1);
    }

    
}

static char * name_unextracted_config_file(const unsigned config_nr, char *unextracted_config_filename) {
    sprintf(unextracted_config_filename, "%s%s_%dx%d/Run1_cfg_%d.lime", configs_dir_name_in, config_template, N_SPC, N_T, config_nr);

    strcat(unextracted_config_filename, extension_config_in);

    return unextracted_config_filename;
}

static void extract_config(const unsigned config_nr, const char *restrict config_filename) {
    char unextracted_config_filename[MAX_LENGTH_NAME];
    char command_lime[MAX_LENGTH_NAME];

    name_unextracted_config_file(config_nr, unextracted_config_filename);

    sprintf(command_lime, "./lime_extract_record %s 2 4 %s", name_unextracted_config_file(config_nr, unextracted_config_filename), config_filename);
    printf("Launching command: %s\n", command_lime);
    printf("This will create a file named %s\n", config_filename);

    if (system(command_lime)) {
        fprintf(stderr, "Problem extracting config %d.\n", config_nr);
        exit(1);
    }
}

static char * name_configuration_file(const unsigned config_nr, char *config_filename) {
    sprintf(config_filename, "%s%s_%dx%d_%d", configs_dir_name_out, config_template, N_SPC, N_T, config_nr);
    
    return config_filename;
}



void SU3_load_config(const unsigned config_nr, mtrx_3x3 *U) {
    bool exit_status = 0;
    //	Loads a link configuration from the file with filename to U.

    in_cfg_data_type *U_in;

#ifdef NEED_CONV_TO_WORKING_PRECISION

    U_in = (in_cfg_data_type *)calloc(VOLUME * DIM, sizeof(in_cfg_data_type));
    TEST_ALLOCATION(U_in);

#else

    U_in = (in_cfg_data_type *)U;

#endif

    char config_filename[MAX_LENGTH_NAME];

    strcat(name_configuration_file(config_nr, config_filename), extension_config_in);
    char mv_command[MAX_LENGTH_NAME+100];
    sprintf(mv_command, "mv %s ./configs/%dx%d/", "./configs/Gen2_24x16_1000.cfg", N_SPC, N_T);
    system(mv_command);

    // extract_config(config_nr, strcat(name_configuration_file(config_nr, config_filename), extension_config_in));

    FILE *config_file;

    printf("Loading: %s.\n", config_filename);
    if ((config_file = fopen(config_filename, "rb")) == NULL) {
        fprintf(stderr, "Error opening config from file %s.\n", config_filename);
        
        exit(EXIT_FAILURE);
    }

    if (fread(U_in, sizeof(in_cfg_data_type), VOLUME * DIM, config_file) == VOLUME * DIM) {
        printf("U Loaded OK for config %d.\n", config_nr);
    } else {
        fprintf(stderr, "Configuration loading failed for config %d.\n", config_nr);
        exit(EXIT_FAILURE);
    }

    fgetc(config_file);

    if (!feof(config_file)) {
        fprintf(stderr, "File has not been read till the end. Check lattice sizes.\n");
        exit(EXIT_FAILURE);
    }

    fclose(config_file);
    sprintf(mv_command, "mv ./configs/%dx%d/%s ./configs/", N_SPC, N_T, "Gen2_24x16_1000.cfg");
    system(mv_command);
    // printf("Removing %s\n", config_filename);
    // remove(config_filename);

#ifdef NEED_BYTE_SWAP_IN

    byte_swap(U_in, sizeof(in_data_type) / 2, VOLUME * DIM * sizeof(in_cfg_data_type));

#endif

#ifdef NEED_CONV_TO_WORKING_PRECISION
    SU3_convert_config_in_work(U_in, U);
    free(U_in);
#endif
}

void SU3_write_config(const unsigned config_nr, mtrx_3x3 *U) {
    //  Loads a link configuration from the file with filename to U.

    out_cfg_data_type *U_out;

#ifdef NEED_CONV_FROM_WORKING_PRECISION

    U_out = (out_cfg_data_type *)calloc(VOLUME * DIM, sizeof(out_cfg_data_type));
    TEST_ALLOCATION(U_out);

    SU3_convert_config_work_out(U, U_out);

#else

    U_out = (out_cfg_data_type *)U;

#endif

#ifdef NEED_BYTE_SWAP_OUT

    byte_swap(U_out, sizeof(out_data_type) / 2, VOLUME * DIM * sizeof(out_cfg_data_type));

#endif

    char config_filename[MAX_LENGTH_NAME];
    name_configuration_file(config_nr, config_filename);

    FILE *config_file;

    printf("Creating: %s.\n", strcat(config_filename, extension_config_out));
    if ((config_file = fopen(config_filename, "wb")) == NULL) {
        fprintf(stderr, "Error creating file %s for config.\n", config_filename);
        exit(EXIT_FAILURE);
    }

    if (fwrite(U_out, sizeof(out_cfg_data_type), VOLUME * DIM, config_file) == VOLUME * DIM) {
        printf("U written OK for config %d.\n", config_nr);

    } else {
        fprintf(stderr, "Configuration writing failed for config %d.\n", config_nr);
    }

    fclose(config_file);

#ifdef NEED_CONV_FROM_DOUBLE
    free(U_out);
#endif
}


static char * name_gaugetransf_file(const unsigned config_nr, char *gaugetransf_filename) {
    sprintf(gaugetransf_filename, "%s%s_%d", configs_dir_name_out, config_template, config_nr);
    
    return gaugetransf_filename;
}


void SU3_write_gauge_transf(const unsigned config_nr, mtrx_3x3 *G) {
    //  Loads a link configuration from the file with filename to U.

    out_cfg_data_type *G_out;

#ifdef NEED_CONV_FROM_WORKING_PRECISION

    G_out = (out_cfg_data_type *)calloc(VOLUME, sizeof(out_cfg_data_type));
    TEST_ALLOCATION(G_out);

    SU3_convert_gaugetransf_work_out(G, G_out);

#else

    G_out = (out_cfg_data_type *)G;

#endif

#ifdef NEED_BYTE_SWAP_OUT

    byte_swap(G_out, sizeof(out_data_type) / 2, VOLUME * sizeof(out_cfg_data_type));

#endif

    char gaugetransf_filename[MAX_LENGTH_NAME];
    name_gaugetransf_file(config_nr, gaugetransf_filename);

    FILE *gaugetransf_file;

    printf("Creating: %s.\n", strcat(gaugetransf_filename, extension_gt_out));
    if ((gaugetransf_file = fopen(gaugetransf_filename, "wb")) == NULL) {
        fprintf(stderr, "Error creating file %s for gauge transformation.\n", gaugetransf_filename);
        exit(EXIT_FAILURE);
    }

    if (fwrite(G_out, sizeof(out_cfg_data_type), VOLUME, gaugetransf_file) == VOLUME ) {
        printf("G written OK for config %d.\n", config_nr);

    } else {
        fprintf(stderr, "Gauging transformation writing failed for config %d.\n", config_nr);
    }

    fclose(gaugetransf_file);

#ifdef NEED_CONV_FROM_DOUBLE
    free(G_out);
#endif
}



void SU3_load_gauge_transf(const unsigned config_nr, mtrx_3x3 * restrict G) {
    bool exit_status = 0;
    //	Loads a link configuration from the file with filename to U.

    in_cfg_data_type *G_in;

#ifdef NEED_CONV_TO_WORKING_PRECISION

    G_in = (in_cfg_data_type *)calloc( VOLUME , sizeof(in_cfg_data_type));
    TEST_ALLOCATION(G_in);

#else

    G_in = (in_cfg_data_type *)G;

#endif

    char gaugetransf_filename[MAX_LENGTH_NAME];

    strcat(name_gaugetransf_file(config_nr, gaugetransf_filename), extension_gt_in);
    
    FILE *gaugetransf_file;

    printf("Loading: %s.\n", gaugetransf_filename);
    if ((gaugetransf_file = fopen(gaugetransf_filename, "rb")) == NULL) {
        fprintf(stderr, "Error opening gauge transformation from file %s.\n", gaugetransf_filename);
        
        exit(EXIT_FAILURE);
    }

    if (fread(G_in, sizeof(in_cfg_data_type), VOLUME, gaugetransf_file) == VOLUME ) {
        printf("G Loaded OK for config %d.\n", config_nr);
    } else {
        fprintf(stderr, "Gauge transformation loading failed for config %d.\n", config_nr);
        exit(EXIT_FAILURE);
    }

    fgetc(gaugetransf_file);

    if (!feof(gaugetransf_file)) {
        fprintf(stderr, "File has not been read till the end. Check lattice sizes.\n");
        exit(EXIT_FAILURE);
    }

    fclose(gaugetransf_file);
    

#ifdef NEED_CONV_TO_WORKING_PRECISION
    SU3_convert_config_in_work(G_in, G);
    free(G_in);
#endif
}

