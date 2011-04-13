/*
 * This is to workaround the buggy sqrt(double) on Windows (VC++ and GCC).
 * This file must be compiled with optimization, as a C file (not C++),
 * then linked into libCGAL: cl -c -O2 -nologo workaround_4_ms.c
 *
 * Sylvain.Pion@sophia.inria.fr, October 1999.
 */

#include <math.h>

double CGAL_ms_sqrt(double d)
{
    return sqrt(d);
}
