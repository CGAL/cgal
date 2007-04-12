// This example program solves the quadratic program
//    minimize    x^T D x + c^T x + c0
//    subject to  A x <rel> b
//                L <= x <= u
// where <rel> stands for a vector of relations from {<=, ==, >=}
//
// In the concrete example below, we have the following data:
//     D          c         A          b      L         U          c0
// -----------------------------------------------------------------------
//     1.0 0.0    0.0 3.0   0.5 1.0    0.5    0.0 0.0   inf, inf    1
//     0.0 0.0
//
// In other words, we minimize x^2 + 3y subject to x/2 + y >= 1/2
// and x,y >= 0. Let's see what this gives: it is clear that in an
// optimal solution, the constraint is satisfied with equality,
// otherwise, we could decrease the objective function value. Let's
// substitute x with 1 - 2y in the objective function; we get
// 4y^2 - y + 1. The derivative is 8y-1, so this is minimzed
// for y = 1/8. This value (as well as x = 1 - 2(1/8) = 3/4)
// are within the bounds, so the solution should be x = 3/4
// and y = 1/8. The objective function value should be
// 9/16 + 3/8 + 1 = 31/16

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#endif

// program and solution types
typedef CGAL::Nonnegative_quadratic_program_from_pointers<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  // program data
  double A_col_0[] = {0.5};
  double A_col_1[] = {1.0};
  double *cols_of_A[] = {A_col_0, A_col_1};
  double b[]       = {0.5};
  CGAL::Comparison_result rt[] = {CGAL::LARGER};
  double D_row_0[] = {2.0, 0.0};
  double D_row_1[] = {0.0, 0.0};
  double *rows_of_D[] = {D_row_0, D_row_1};
  double c[]       = {0.0, 3.0};
  double c0 = 1.0;

  // now construct the quadratic program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program qp (2, 1, cols_of_A, b, rt, rows_of_D, c, c0);
  // solve the program
  Solution s = CGAL::solve_nonnegative_quadratic_program(qp, ET());

  // output solution
  if (s.status() == CGAL::QP_OPTIMAL) { // we know that, don't we?
    std::cout << "Optimal feasible solution: ";
    for (Solution::Variable_value_iterator it = s.variable_values_begin();
	 it != s.variable_values_end(); ++it)
      std::cout << CGAL::to_double(*it) << "  ";
    std::cout << "\nOptimal value: " <<  CGAL::to_double(s.objective_value())
	      << std::endl;
  }

  return 0;
}
