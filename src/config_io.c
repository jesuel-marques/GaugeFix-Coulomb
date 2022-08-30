#include <complex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <string.h>

#include <SU3_gaugefixing_parameters.h>  //	Gauge-fixing specific parameters
#include <SU3_parameters.h>              //	Simulation parameters
#include <SU2_ops.h>                        //	SU(2) operations
#include <SU3_ops.h>                        //	SU(3) operations
                                            //	positions and links on the lattice.

#include <lattice.h>  //	Initialization functions and calculations of

char configs_dir_name_in[200];	//	input from command line
char configs_dir_name_out[200];	//	input from command line
char config_template[10] ;	//	input from command line

char gaugetransf_dir_name_in[200];
char gaugetransf_dir_name_out[200];

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

void greeter_function(const char * restrict program_name) {
    printf("Hello %s!\n.", getenv("USER"));
    printf("Program %s compiled at %s on %s\n", program_name, __TIME__, __DATE__);
    printf("Using C version: %ld\n", __STDC_VERSION__);
    // printf
    // #include "../README.txt"
    // ;
    printf("\n");

}

void handle_input(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stdout, "Usage: Input names for input and output config directory and basename.\n");

        if (argc > 4) {
            fprintf(stdout, " only");
        }
    }

    if (argc != 4) {
        fprintf(stderr, ".\n");
        exit(EXIT_FAILURE);
    }

    char size_directory[20];

    sprintf(size_directory, "%dx%d/", N_SPC, N_T);

    strcpy(configs_dir_name_in,  argv[1]);
    strcpy(configs_dir_name_out, argv[2]);

    strcpy(gaugetransf_dir_name_in , argv[1]);
    strcpy(gaugetransf_dir_name_out, argv[2]);

    strcat(configs_dir_name_in , size_directory);
    strcat(configs_dir_name_out, size_directory);

    strcat(gaugetransf_dir_name_in , size_directory);
    strcat(gaugetransf_dir_name_out, size_directory);

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

int create_output_directory(void) {
    printf("Creating %s.\n", configs_dir_name_out);
    
    char command[MAX_LENGTH_NAME * 2];
    sprintf(command, "test -d %s", configs_dir_name_out);
    
    switch (system(command)){
        case 0:
            printf("Warning: Directory already exists.\n");
            return 0;

        case 1:
            sprintf(command, "mkdir %s", configs_dir_name_out);
            int exit_status = system(command);
            
            if(exit_status){
                fprintf(stderr, "Error: Some problem occured when creating output directory.\n");
                fprintf(stderr, "Exit code of command '%s': %d.\n", command, exit_status);
            }
            return exit_status;

            
        default:
            fprintf(stderr, "Error: Problem with command %s. Could not test existence of directory.\n", command);
            return -1;
    }

    return -1;
}

static char * name_unextracted_config_file(const unsigned config_nr, char * restrict unextracted_config_filename) {
    int exit_status = sprintf(unextracted_config_filename, "%sRun1_cfg_%d.lime", configs_dir_name_in, config_nr);
    
    if(exit_status < 0){
        return '\0';
    }
    else{
        printf("Unextracted config file name is %s\n", unextracted_config_filename);
        return unextracted_config_filename;
    }
}

static int extract_config(const unsigned config_nr, const char * restrict config_filename) {
    char unextracted_config_filename[MAX_LENGTH_NAME];
    char command_lime[MAX_LENGTH_NAME+1000];

    name_unextracted_config_file(config_nr, unextracted_config_filename);
    

    int exit_status = sprintf(command_lime, "/ichec/home/users/jesuel/chroma-install/install-files/qdp++/bin/lime_extract_record %s 2 4 %s", 
                                name_unextracted_config_file(config_nr, unextracted_config_filename), config_filename);

    // int exit_status = sprintf(command_lime, "./lime_extract_record %s 2 4 %s", 
    //                             name_unextracted_config_file(config_nr, unextracted_config_filename), config_filename);

    if(exit_status < 0){
        fprintf(stderr, "Error: Problem creating lime command.\n");
        return -1;
    }
    
    printf("Launching command: %s\n", command_lime);
    printf("This will create a file named %s\n", config_filename);
    
    if (system(command_lime)) {

        return -1;

    }
    
    return 0;
    
}

static char * name_configuration_file(const unsigned config_nr, char * restrict config_filename) {
    int exit_status = sprintf(config_filename, "%s%s_%dx%d_%d", configs_dir_name_out, config_template, N_SPC, N_T, config_nr);
    if(exit_status < 0){
        return '\0';
    }
    else{    
        return config_filename;
    }
}



int SU3_load_config(const unsigned config_nr, mtrx_3x3 * restrict U) {
    //	Loads a link configuration from the file with filename to U.

    char config_filename[MAX_LENGTH_NAME];

    if( name_configuration_file(config_nr, config_filename) == '\0'){
        fprintf(stderr, "Error: Problem naming configuration.\n");
        return -1;
    }

    strcat(config_filename, extension_config_in);

    // char mv_command[MAX_LENGTH_NAME*3];
    // sprintf(mv_command, "mv %s ./configs/%dx%d/", "./configs/Gen2_24x16_1000.cfg", N_SPC, N_T);
    // system(mv_command);

    if(extract_config(config_nr, strcat(name_configuration_file(config_nr, config_filename), extension_config_in))){
        fprintf(stderr, "Error: Problem extracting config %d.\n", config_nr);
        return -1;
    }
    printf("Config %d extracted OK.\n", config_nr);

    FILE *config_file;

    if ((config_file = fopen(config_filename, "rb")) == NULL) {

        fprintf(stderr, "Error: Problem opening config file %s.\n", config_filename);

        return -1;
    }

    in_cfg_data_type *U_in;

    #ifdef CONV_CFG_TO_WORKING_PRECISION

        U_in = (in_cfg_data_type *)calloc(VOLUME * DIM, sizeof(in_cfg_data_type));
        if(TEST_ALLOCATION(U_in)){

            return -1;
        }

    #else

        U_in = (in_cfg_data_type *)U;

    #endif
    
    printf("Loading: %s.\n", config_filename);

    if (fread(U_in, sizeof(in_cfg_data_type), VOLUME * DIM, config_file) != VOLUME * DIM) {

        fprintf(stderr, "Error: Reading from file failed for config %d.\n", config_nr);
        #ifdef CONV_CFG_TO_WORKING_PRECISION 
             free(U_in);
        #endif
        fclose(config_file);
        return -1;

    }

    fgetc(config_file);

    if (!feof(config_file)) {
        fprintf(stderr, "Error: File has not been read till the end. Check lattice sizes.\n");
        #ifdef CONV_CFG_TO_WORKING_PRECISION 
             free(U_in);
        #endif
        fclose(config_file);
        return -1;
    }

    fclose(config_file);
    // sprintf(mv_command, "mv ./configs/%dx%d/%s ./configs/", N_SPC, N_T, "Gen2_24x16_1000.cfg");
    // system(mv_command);
    printf("Removing %s\n", config_filename); 
    if(remove(config_filename)){
        fprintf(stderr, "Error: Problem removing %s.\n", config_filename);
    }
    else{
        printf("%s succesfully removed.\n", config_filename);
    }

    #ifdef NEED_BYTE_SWAP_IN

        if(byte_swap(U_in, sizeof(in_data_type) / 2, VOLUME * DIM * sizeof(in_cfg_data_type))){
            fprintf(stderr, "Error: Problem with the byte_swap. Unknown size.\n");
            #ifdef CONV_CFG_TO_WORKING_PRECISION
                free(U_in);
            #endif
            return -1;
        }

    #endif

    #ifdef CONV_CFG_TO_WORKING_PRECISION
        SU3_convert_cfg_in_work(U_in, U);
        free(U_in);
    #endif

    return 0;
}

int SU3_write_config(const unsigned config_nr, mtrx_3x3 * restrict U) {
    //  Loads a link configuration from the file with filename to U.

    out_cfg_data_type *U_out;

    #ifdef CONV_CFG_FROM_WORKING_PRECISION

        U_out = (out_cfg_data_type *)calloc(VOLUME * DIM, sizeof(out_cfg_data_type));
        TEST_ALLOCATION(U_out);

        SU3_convert_cfg_work_out(U, U_out);

    #else

        U_out = (out_cfg_data_type *)U;

    #endif

    #ifdef NEED_BYTE_SWAP_OUT
        if (byte_swap(U_out, sizeof(out_data_type) / 2, VOLUME * DIM * sizeof(out_cfg_data_type))){
            fprintf(stderr, "Error: Problem with the byte_swap. Unknown size.\n");
            #ifdef CONV_CFG_FROM_WORKING_PRECISION
                free(U_out);
            #endif
            return -1;
        }

    #endif

    char config_filename[MAX_LENGTH_NAME];
    name_configuration_file(config_nr, config_filename);

    FILE *config_file;

    printf("Creating: %s.\n", strcat(config_filename, extension_config_out));
    if ((config_file = fopen(config_filename, "wb")) == NULL) {

        fprintf(stderr, "Error: Problem creating file %s for config.\n", config_filename);

        #ifdef CONV_CFG_FROM_WORKING_PRECISION
            free(U_out);
        #endif

        return -1;

    }

    if (fwrite(U_out, sizeof(out_cfg_data_type), VOLUME * DIM, config_file) != VOLUME * DIM) {

        fclose(config_file);

        #ifdef CONV_CFG_FROM_WORKING_PRECISION
            free(U_out);
        #endif

        return -1;
    }

    fclose(config_file);

    #ifdef CONV_CFG_FROM_WORKING_PRECISION
        free(U_out);
    #endif

    return 0;
}


static char * name_gaugetransf_file(const unsigned config_nr, char * restrict gaugetransf_filename) {
    int exit_status = sprintf(gaugetransf_filename, "%s%s_%d", gaugetransf_dir_name_out, config_template, config_nr);
    
    if(exit_status < 0){

        return '\0';

    }
    else{

        return gaugetransf_filename;

    }
}


int SU3_write_gauge_transf(const unsigned config_nr, mtrx_3x3 * restrict G) {
    //  Loads a link configuration from the file with filename to U.

    out_gt_data_type *G_out;

    #ifdef CONV_GT_FROM_WORKING_PRECISION

        G_out = (out_gt_data_type *)calloc(VOLUME, sizeof(out_gt_data_type));
        if(TEST_ALLOCATION(G_out)){
            
            return -1;
        
        }
        SU3_convert_gt_work_out(G, G_out);

    #else

        G_out = (out_gt_data_type *)G;

    #endif

    char gaugetransf_filename[MAX_LENGTH_NAME];
    name_gaugetransf_file(config_nr, gaugetransf_filename);

    FILE *gaugetransf_file;

    printf("Creating: %s.\n", strcat(gaugetransf_filename, extension_gt_out));
    if ((gaugetransf_file = fopen(gaugetransf_filename, "wb")) == NULL) {

        fprintf(stderr, "Error: Problem opening file %s to store gauge transformation.\n", gaugetransf_filename);
        
        #ifdef CONV_GT_FROM_WORKING_PRECISION

            free(G_out);
        
        #endif

        return -1;
    }

    if (fwrite(G_out, sizeof(out_gt_data_type), VOLUME, gaugetransf_file) != VOLUME ) {

        fprintf(stderr, "Error: Problem writing gauge transformation to file %s.\n", gaugetransf_filename);

        fclose(gaugetransf_file);

        #ifdef CONV_GT_FROM_WORKING_PRECISION

            free(G_out);
            
        #endif

        return -1;
    }

    fclose(gaugetransf_file);

    #ifdef CONV_GT_FROM_WORKING_PRECISION

        free(G_out);
        
    #endif

    return 0;
}

int SU3_load_gauge_transf(const unsigned config_nr, mtrx_3x3 * restrict G) {
    //	Loads a gauge transformation to G.

    in_gt_data_type *G_in;

    #ifdef CONV_GT_TO_WORKING_PRECISION

        G_in = (in_gt_data_type *)calloc( VOLUME , sizeof(in_gt_data_type));
        if(TEST_ALLOCATION(G_in)){
            return -1;
        }

    #else

        G_in = (in_gt_data_type *)G;

    #endif

    char gaugetransf_filename[MAX_LENGTH_NAME];

    if(name_gaugetransf_file(config_nr, gaugetransf_filename) == '\0'){
        fprintf(stderr, "Error: Problem naming gauge transformation.\n");
        
        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif

        return -1;
    }


    strcat(gaugetransf_filename, extension_gt_in);
    
    FILE *gaugetransf_file;

    printf("Loading: %s.\n", gaugetransf_filename);
    if ((gaugetransf_file = fopen(gaugetransf_filename, "rb")) == NULL) {

        fprintf(stderr, "Error: Problem opening file %s to load gauge transformation.\n", gaugetransf_filename);

        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif

        return -1;

    }

    if (fread(G_in, sizeof(in_gt_data_type), VOLUME, gaugetransf_file) != VOLUME ) {

        fprintf(stderr, "Error: Problem reading file %s to load gauge transformation.\n", gaugetransf_filename);

        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif

        fclose(gaugetransf_file);
        return -1;

    }

    fgetc(gaugetransf_file);

    if (!feof(gaugetransf_file)) {

        fprintf(stderr, "Error: File has not been read till the end. Check lattice sizes.\n");
        
        #ifdef CONV_GT_TO_WORKING_PRECISION
            free(G_in);
        #endif

        fclose(gaugetransf_file);
        return -1;
    }

    fclose(gaugetransf_file);

    #ifdef CONV_GT_TO_WORKING_PRECISION
        SU3_convert_gt_in_work(G_in, G);
        free(G_in);
    #endif

    return 0;
}