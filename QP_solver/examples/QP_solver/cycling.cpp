// example: solve a linear program that by default leads to cycling,
// using Bland pricing
#include <iostream>
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// choose exact floating-point type
#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#endif

// program and solution types
typedef CGAL::Quadratic_program_from_mps<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  std::ifstream in ("cycling.mps");
  Program lp(in);         // read program from file
  assert (lp.is_valid()); // we should have a valid mps file

  // solve the program, using ET as the exact type
  // choose verbose mode and Bland pricing
  CGAL::Quadratic_program_options options;
  options.set_verbosity(1);                              // verbose mode 
  options.set_pricing_strategy(CGAL::QP_EXACT_BLAND);    // Bland's rule
  Solution s = CGAL::solve_nonnegative_linear_program(lp, ET(), options);

  // output solution
  if (s.is_optimal()) { // we know that, don't we?
    std::cout << "Optimal feasible solution: ";
    for (Solution::Variable_value_iterator it = s.variable_values_begin();
	 it != s.variable_values_end(); ++it)
      std::cout << *it << "  ";
    std::cout << std::endl << "Optimal objective function value: "
	      << s.objective_value() << std::endl;
  }

  return 0;
}
