// example: construct a nonnegative linear program from data
// the LP below is the first nonnegative linear program example
// in the user manual
#include <iostream>
#include <cassert>

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
typedef CGAL::Quadratic_program<int> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  // by default, we have a nonnegative LP with Ax <= b
  Program lp (CGAL::SMALLER, true, 0, false, 0);

  // now set the non-default entries
  const int X = 0;
  const int Y = 1;
  lp.set_a(X, 0,  1); lp.set_a(Y, 0, 1); lp.set_b(0, 7);  //  x + y  <= 7
  lp.set_a(X, 1, -1); lp.set_a(Y, 1, 2); lp.set_b(1, 4);  // -x + 2y <= 4
  lp.set_c(Y, -32);                                       // -32y

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_linear_program(lp, ET());
  assert (s.solves_nonnegative_linear_program(lp));

  // output solution
  std::cout << s;
  return 0;
}
