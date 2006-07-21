#ifndef CGAL_GBRS_SOLVE_1_H
#define CGAL_GBRS_SOLVE_1_H

#include <mpfi.h>
#include <CGAL/Gbrs_polynomial_1.h>

CGAL_BEGIN_NAMESPACE
// both functions return the number of roots, -1 if it is infinite and -2 if
// there were an error

// solve given the preciseness
int solve_1 (mpfi_t *&, const Rational_polynomial_1 &, unsigned int);

// solve with the default preciseness
int solve_1 (mpfi_t *&, const Rational_polynomial_1 &);

CGAL_END_NAMESPACE

#endif	// CGAL_GBRS_SOLVE_1_H
