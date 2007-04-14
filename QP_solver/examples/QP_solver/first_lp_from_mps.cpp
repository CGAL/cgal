#include <iostream>
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// choose exact integral type
#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#endif

// program and solution types
typedef CGAL::Quadratic_program_from_mps<int> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  std::ifstream in ("first_lp.mps");
  Program lp(in);         // read program from file
  assert (lp.is_valid()); // we should have a valid mps file,...
  assert (lp.is_linear());// ...and it encodes a linear program

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_linear_program(lp, ET());

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
