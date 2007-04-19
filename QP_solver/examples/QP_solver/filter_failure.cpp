// example: solve a nonnegative quadratic program that by default leads 
// to double overflows, using Dantzig pricing
#include <iostream>
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// choose exact floating-point type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif
// program and solution types
typedef CGAL::Quadratic_program_from_mps<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  std::ifstream in ("filter_failure.mps");
  Program qp(in);         // read program from file
  assert (qp.is_valid()); // we should have a valid mps file
  assert (qp.is_nonnegative()); // .. and it should be nonnegative

  // solve the program, using ET as the exact type
  // choose verbose mode and non-filtered Dantzig pricing
  CGAL::Quadratic_program_options options;
  options.set_pricing_strategy(CGAL::QP_DANTZIG);        // Dantzig's rule
  Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET(), options);

  // output solution
  std::cout << s;
  return 0;
}
