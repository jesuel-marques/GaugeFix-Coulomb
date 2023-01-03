#ifndef MISC_H
#define MISC_H

#define OMP_PARALLEL_BASIC "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
#define OMP_PARALLEL_FOR  _Pragma(OMP_PARALLEL_BASIC)

#define POW2(a) (a) * (a)   //  squares a number

#endif  //MISC_H