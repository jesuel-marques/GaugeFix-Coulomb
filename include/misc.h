/*
    Simple miscellaneous operations.

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

#ifndef MISC_H
#define MISC_H

/* Used to shorten the pragma directive which paralellizes loops */
#define OMP_PARALLEL_BASIC "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
/* Used to shorten the pragma directive which paralellizes loops */
#define OMP_PARALLEL_FOR  _Pragma(OMP_PARALLEL_BASIC)

/* squares a number */
#define POW2(a) (a) * (a)

#endif  //MISC_H