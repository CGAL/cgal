// example: construct a nonnegative quadratic program from given iterators
// the QP below is the first nonnegative quadratic program example
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
typedef CGAL::Nonnegative_quadratic_program_from_iterators
<int**,                                                // for A
 int*,                                                 // for b
 CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
 int**,                                                // for D
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
  int  D1[] = {2};                            // 2D_{1,1}
  int  D2[] = {0, 8};                         // 2D_{2,1}, 2D_{2,2}
  int*  D[] = {D1, D2};                       // D-entries on/below diagonal
  int   c[] = {0, -32};
  int  c0   = 64;                             // constant term

  // now construct the quadratic program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program qp (2, 2, A, b, r, D, c, c0);

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET());

  // output solution
  std::cout << s;
  return 0;
}
