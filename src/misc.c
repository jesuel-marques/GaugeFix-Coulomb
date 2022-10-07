#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <misc.h>
#include <SU3_parameters.h>

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

void greeter_function(const char * restrict program_name) {
    printf("Hello %s!\n", getenv("USER"));
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

short test_allocation_function(const void *pointer, const char *location) {
    //	Test if allocation was successful.
    if (pointer == NULL) {
        fprintf(stderr, "Memory allocation failed at %s.\n", location);
        return -1;
    }

    return 0;
}