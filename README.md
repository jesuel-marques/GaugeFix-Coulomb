# GaugeFix-Coulomb

This program gauge-fixes lattice QCD gauge configurations to Coulomb gauge.

## Install

### On Linux:

Go to the main directory, where gauge_fix_coulomb.c is located. 
The command to compile is 
```
make
```
Make will read the Makefile and should know how to compile the code. The Makefile 
provided is set to use GCC as the C compiler. If you use something else, you would
have to change that and change the flags appropriately. Also feel free to change the 
flags in the Makefile to suit your needs.

### Other platforms:

Not guaranteed to work, but a C compiler should be able to compile the code anyway.

## Configuration

Before you run the program, you need know if:

- Your configuration is in the SCIDAC format and if the xml header has been stripped 
away: 
This is currently the only supported format. This is the standard output format from 
USQCD-chroma. To strip away the xml header, use the command
```
lime_extract_record path_to_configs/config.lime 2 4 path_to_extracted_configs/config.cfg
```
If you have chroma installed, you should be able to run this command.

It will leave the binary data which the program reads, and which contains only the 
configuration itself, and no metadata. The size of this file should be
    nt * ns * ns * ns * dimensions * Nc * Nc * 2 * sizeof(precision) bytes

where dimensions = 4, Nc = 3 and sizeof(precision) is 4 for single (float) and 8 for 
double precision, ns is the number of points in the spatial directions and nt is the 
number of points in the temporal direction. If your file is longer or shorter than this, 
there is something wrong.

- If this is not the format you have, and you still want to use this program, you must 
use another program to convert the configurations to the supported format first or write
a module that does that internally. 

- Your configuration is little or big endian: If you need to change the endianness, 
you can do so by leaving the line 
```
#define NEED_BYTE_SWAP_IN
```
uncommented in include/flags.h. You can also change the endianess of an output 
configuration if you leave
```
#define NEED_BYTE_SWAP_OUT
```
uncommented in include/flags.h. 

You cannot change the endianness of the gauge-transformation, neither for input or 
output.

- Your configurations are in single or double precision:

All internal calculations are in double precision. If your input configurations are in 
single precision, leave the line 
```
#define CONV_CFG_TO_WORKING_PRECISION
```
uncommented in include/flags.h. 

If you want output configurations to be written in float, leave 
```
#define CONV_CFG_FROM_WORKING_PRECISION
```
uncommented in include/flags.h. 

If you want to input or output gauge-fixing transformations, also leave the following 
lines commented or uncommented as appropriate:
```
#define  CONV_GT_TO_WORKING_PRECISION
#define  CONV_GT_FROM_WORKING_PRECISION
```
following the same pattern as for the configurations.

## Usage

The main executable is gauge_fix_coulomb. It is a program which reads in a gauge
configuration in the SCIDAC binary format, together with some other parameters and 
outputs a gauge transformation which brings the configuration to Coulomb-gauge.

In the main directory, run

./gauge_fix_coulomb config_path gtransformation_path ns nt parameters.par

where

config_path: the path to your configuration in need of gauge-fixing, which is an input;
gtransformation_path: the path to a gauge-transformation which brings your configuration
to Coulomb-gauge, which is the main output;
ns: number of points in the spatial directions;
nt: number of points in the temporal direction;
parameters.par: parameter file, which sets parameters for the gauge-fixing.

On a successful execution, it should output the gauge transformation to 
gtransformation_path and also a file "sweeps_to_gaugefix_nsxnt.txt", which records
how many sweeps were needed to gauge-fix the particular configuration. It writes the
gtransformation_path in the first column and the number of sweeps in the next column.

## Parameter file format

You can specify a few parameters for the gauge-fixing. If a parameter is not specified,
the program will use some default value. The parameters are specified in a simple text
file with one parameter per line in a "key: value" format.

The algorithm used for the gauge-fixing does sweeps on the configurations and at each
sweep calculates an update for the gauge-fixing field. This update is such that you are
always getting closer to the gauge-fixing condition.

- In order to monitor how close to the gauge-fixing condition one is, we use the field
gfix_proxy. This calculates a scalar quantity which is 0 if the gauge-fixing condition 
is respected. Currently there are two such gfix_proxy implemented, which are called e2 
and theta.

theta = 1/(N_c V) sum_n Re Tr[(div.A(n)).(div.A(n))^dagger]

and e2 is the normalized sum of the squares of the color components of the divergence of
A, in which the spatial divergence is meant. 
The default value for gfix_proxy is e2.

- You can set how close you want to get with the 'tolerance' key. 
The default value is 1E-16. It has to be a positive number.

- You can set a maximum number of sweeps that you are willing to wait for the program 
to gauge-fix the configuration with the 'max_sweeps_to_fix' key. 
The default value is 100000. It has to be a positive number.

- The program doesn't check if you arrived at the specified tolerance all the time, 
because this is wasteful. You can provide an estimate of when it should first check with
the 'estimate_sweeps_to_gfixprogress' key. This is roughly the number of sweeps it will
wait to check if you got to the specified tolerance. After that, the program will check
more frequently, as you get closer and closer. 
The default value is 1000. It has to be a positive number.

- It is important to keep the condition that the links are SU(3) elements. In order to do 
this, the program applies a reunitarization step to correct for accumulated numerical 
errors. You can set the number of sweeps to wait until you do a reunitarization step 
with the 'sweeps_to_reunitarization'.
The default value is 250. It has to be a positive number.

- There are two parameters that are specific to the overrelaxation algorithm. Their 
meaning will be explained in the in the Internal workings section. 

max_hits: Should be an integer.
omega_OR: Should be a number 1 <= omega_OR < 2.

Feel free to play around with these values, but in my experience, putting max_hits to 2
and omega_OR to a number close to 2, like 1.97 is the fastest. 

If you put max_hits to anything above 2 and keep omega_OR the same, it just takes longer
to do the same job. If you put max_hits to 1, it takes much longer, and the same happens
if you put omega_OR to a lower value. If you see a different behaviour, let me know. The
default value for max_hits is 2 and for omega_OR is 1.95.

You can leave blank lines between parameters and you can have as many whitespaces as you
want between the key the ':' delimiter and the value. 

If you supply a parameter which is not supported, the program will warn about this and 
will ignore the parameter. 

It will also ignore any line that does not have the "key:value" format. So you can use 
this fact to make comments in the parameter file. However, if your comment contains ":",
it will think that you wanted somehow to specify a parameter with a weird name and warn 
you, but that does no harm.

## Internal workings

In order to fix the gauge the overrelaxation algorithm. To know more about this take a 
look at DOI: 10.22323/1.396.0057.

From a user point of view, you need to specify two parameters: max_hits and omega_OR.
These parameters have to do with two steps of the algorithm. The gauge-fixing condition
translates on the lattice to a maximization of a global functional. However, as we don't
know how to maximize this global functional, we must resort to a local maximization. 

The first step of the algorithm has to do with this local maximization. At each lattice 
site, we should find g, an SU(3) matrix, which maximizes

Re Tr [g.w].

For SU(3) this is done via an interactive procedure. The max_hits is the numbers of 
iterations in this process. If you want to get closer to the local maximum, you need to 
increase the number of hits. This does not translate directly into a faster gauge-fixing
however, because the maximum which we actually want is the global maximum and maximizing
everything locally is not the same as the global maximization. Also, these hits take 
some time to perform, so what you may gain from getting closer to the local maximum may 
not compensate for this extra time.

The second step in the algorithm is the actual overrelaxation. Instead of using the g 
that you found in the previous step, you use g^omega_OR, a power of the matrix between 1
and 2. This speeds up the conversion considerably, and makes the scaling of the number
of sweeps with the lattice size increase linearly, instead of quadratically, at least if
you use an omega_OR which is close to the (in principle) unknown optimal value. 
From my experience, for large enough lattices, the optimal value is actually 
surprisingly close to 2.

## Exit codes

0 if successful and a negative integer if some problem occured.

## More on the format

The mapping between the directions and the Lorentz indices is defined in lattice.h.

```
#define X_INDX 0
#define Y_INDX 1
#define Z_INDX 2
#define T_INDX 3
```

So mu=0 corresponds to the spatial x direction, and so on. The indices for transversing
the configurations, are expected to have T as the most external index, and then 
Z, Y, X, Lorentz index, color index, color index and then the real and imaginary parts 
of the matrix. The configuration is read in this order. The full matrices are expected 
to be stored, and not just two rows.

If you need to change the mapping between mu and the directions, you currently need to 
change the define directives. I realize this is not the most convenient and should be 
fixed in the future.

## Not (at least currently) supported

In the order of priority of possible future implementation 

1. [ ] Support for other gauges
2. [ ] Other algorithms for gauge-fixing other than overrelaxation
3. [ ] Other gfix_proxy apart from e2 and theta
4. [ ] Gently pause or stop the running program
5. [ ] Signal handling. Program will terminate abruptly if termination signals are passed
6. [ ] Support for GPU acceleration
7. [ ] Other formats for configurations different from the headerless SCIDAC
8. [ ] Nc different from 3
9. [ ] Lattices with different number of points in different spatial directions
10. [ ] Gauge field in a representation different from the fundamental
11. [ ] Number of dimensions different from 4
12. [ ] Support for custom precision datatypes

Error handling is not great right now, but I intend to improve this as time allows.