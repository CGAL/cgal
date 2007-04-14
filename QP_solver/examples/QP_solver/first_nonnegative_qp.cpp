// example: construct a nonnegative quadratic program from data
// the QP below is the first nonnegative quadratic program example 
// in the user manual
#include <iostream>
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
typedef CGAL::Quadratic_program<int> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  // by default, we have a nonnegative QP with Ax <= b
  Program qp (CGAL::SMALLER, true, 0, false, 0); 
  
  // now set the non-default entries: 0 <-> x, 1 <-> y
  qp.set_a(0, 0,  1); qp.set_a(1, 0, 1); qp.set_b(0, 7);  //  x + y  <= 7
  qp.set_a(0, 1, -1); qp.set_a(1, 1, 2); qp.set_b(1, 4);  // -x + 2y <= 4
  qp.set_d(0, 0, 2); qp.set_d (1, 1, 8);                  // x^2 + 4 y^2
  qp.set_c(1, -32);                                       // -32y
  qp.set_c0(64);                                          // +64

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET());

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
