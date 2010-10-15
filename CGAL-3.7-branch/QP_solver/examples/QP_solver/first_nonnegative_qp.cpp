// example: construct a nonnegative quadratic program from data
// the QP below is the first nonnegative quadratic program example 
// in the user manual
#include <iostream>
#include <cassert>
#include <CGAL/basic.h>
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
  // by default, we have a nonnegative QP with Ax <= b
  Program qp (CGAL::SMALLER, true, 0, false, 0); 
  
  // now set the non-default entries
  const int X = 0; 
  const int Y = 1;
  qp.set_a(X, 0,  1); qp.set_a(Y, 0, 1); qp.set_b(0, 7);  //  x + y  <= 7
  qp.set_a(X, 1, -1); qp.set_a(Y, 1, 2); qp.set_b(1, 4);  // -x + 2y <= 4
  qp.set_d(X, X, 2); qp.set_d (Y, Y, 8); // !!specify 2D!!    x^2 + 4 y^2
  qp.set_c(Y, -32);                                       // -32y
  qp.set_c0(64);                                          // +64

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET());
  assert (s.solves_nonnegative_quadratic_program(qp));
 
  // output solution
  std::cout << s; 
  return 0;
}
