// example: solve a linear program that by default leads to cycling,
// using Bland pricing
#include <iostream>
#include <fstream>

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
  std::ifstream in ("cycling.mps");
  Program lp(in);         // read program from file
  assert (lp.is_valid()); // we should have a valid mps file...
  assert (lp.is_linear());      // ... and it should be linear...
  assert (lp.is_nonnegative()); // as well as nonnegative

  // solve the program, using ET as the exact type
  // choose verbose mode and Bland pricing
  CGAL::Quadratic_program_options options;
  options.set_verbosity(1);                         // verbose mode
  options.set_pricing_strategy(CGAL::QP_BLAND);     // Bland's rule
  options.set_auto_validation(true);                // automatic self-check
  Solution s = CGAL::solve_nonnegative_linear_program(lp, ET(), options);
  assert (s.is_valid());                     // did the self-check succeed?

  // output solution
  std::cout << s;
  return 0;
}
