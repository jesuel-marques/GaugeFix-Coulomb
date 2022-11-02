#include <stdio.h>

#include <misc.h>
#include <SU3_parameters.h>


extern short n_SPC;
extern short n_T;


void greeter_function(const char * restrict program_name) {  
    printf("Hello %s!\n", getenv("USER"));
    printf("Program %s compiled at %s on %s\n", program_name, __TIME__, __DATE__);
    printf("Using C version: %ld\n", __STDC_VERSION__);
    
    printf("\n");

}


short test_allocation_function(const void *pointer, 
                               const char *location) {
    //	Test if allocation was successful.
    if (pointer == NULL) {
        fprintf(stderr, "Memory allocation failed at %s.\n", location);
        return -1;
    }

    return 0;
}


int write_sweeps_to_gaugefix(char * config_filename, 
                             int sweeps){
    FILE* sweeps_to_gaugefix;
    char filename_sweeps_to_gaugefix[MAX_LENGTH_NAME];
    sprintf(filename_sweeps_to_gaugefix, "sweeps_to_gaugefix_%dx%d.txt", n_SPC, n_T);

    if((sweeps_to_gaugefix = fopen(filename_sweeps_to_gaugefix, "a+")) == NULL){

        fprintf(stderr, "Error opening file %s to record sweeps needed to gaugefix.\n", 
                        filename_sweeps_to_gaugefix);
        return -1;

    }
    else{

        fprintf(sweeps_to_gaugefix, "%s\t%d\n", config_filename, sweeps);
        fclose(sweeps_to_gaugefix);

    }

    return 0;

}