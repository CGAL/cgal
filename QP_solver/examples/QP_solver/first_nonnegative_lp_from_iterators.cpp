// example: construct a nonnegative linear program from given iterators
// the LP below is the first nonnegative linear program example
// in the user manual
#include <iostream>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

// program and solution types
typedef CGAL::Nonnegative_linear_program_from_iterators
<int**,                                                // for A
 int*,                                                 // for b
 CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
 int*>                                                 // for c
Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  int  Ax[] = {1, -1};                        // column for x
  int  Ay[] = {1,  2};                        // column for y
  int*  A[] = {Ax, Ay};                       // A comes columnwise
  int   b[] = {7, 4};                         // right-hand side
  CGAL::Const_oneset_iterator<CGAL::Comparison_result>
        r(    CGAL::SMALLER);                 // constraints are "<="
  int   c[] = {0, -32};

  // now construct the linear program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program lp (2, 2, A, b, r, c); // constant term defaults to 0

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_linear_program(lp, ET());

  // output solution
  std:: cout << s;
  return 0;
}
