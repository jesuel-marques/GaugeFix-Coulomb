//	gcc -o plaquette_check plaquette_check.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  source/measurement.c source/config_io.c -lm -O4 -march=skylake -fopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <tgmath.h>


#include <SU3_parameters.h>			//	Simulation parameters

#include <lattice.h>			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include <SU3_ops.h>

#include <matrix_power.h>

const char extension_config_in[]  = ".cfg";
const char extension_config_out[] = "_clmb.cfg";

const char extension_gt_in[] = "_clmb.gt";
const char extension_gt_out[] = "_clmb.gt";

int config_exception_list[] = {-1};

int main(void) {
	mtrx_3x3 a_to_power;
	mtrx_3x3 a = {.m={0.254985810183403  + I * (-0.169733480381236), 0.437294090898545  + I * (-0.638466346569361),-0.508100811963174  + I * ( 0.221677580031758),
					  0.293989643903225  + I * ( 0.393256005929310), 0.188221457593626  + I * ( 0.333447337901830), 0.097870630589627  + I * ( 0.776354750416125),
					  0.708027903935633  + I * (-0.404708626754461),-0.393281432336049  + I * ( 0.315990305374718),-0.275049681392046  + I * (-0.068810822797471)}};

	print_matrix_3x3(&a, "matrix", 16);
	matrix_power_3x3(&a, 0.5, &a_to_power);
	print_matrix_3x3(&a_to_power, "square root of matrix ", 16);
	
	print_matrix_3x3(&a, "matrix", 16);
	matrix_power_3x3(&a, 0.5, &a_to_power);
	print_matrix_3x3(&a_to_power, "square root of matrix ", 16);

	print_matrix_3x3(&a, "matrix", 16);
	matrix_power_3x3(&a, 1/3., &a_to_power);
	print_matrix_3x3(&a_to_power, "cubic root of matrix ", 16);

	return EXIT_SUCCESS;
}