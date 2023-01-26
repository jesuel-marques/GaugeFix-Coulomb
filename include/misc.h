#ifndef MISC_H
#define MISC_H

/* Used to shorten the pragma directive which paralellizes loops */
#define OMP_PARALLEL_BASIC "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
/* Used to shorten the pragma directive which paralellizes loops */
#define OMP_PARALLEL_FOR  _Pragma(OMP_PARALLEL_BASIC)

/* squares a number */
#define POW2(a) (a) * (a)

#endif  //MISC_H