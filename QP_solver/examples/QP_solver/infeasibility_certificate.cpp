// example: extracting and verifying a proof of infeasibility from the solution
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
typedef CGAL::Nonnegative_linear_program_from_iterators
<int**,                                                // for A
 int*,                                                 // for b
 CGAL::Comparison_result*,                             // for r
 int*>                                                 // for c
Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

// we demonstrate Farkas Lemma: either the system
//     A x <= b
//       x >= 0
// has a solution, or there exists y such that
//       y >= 0
//    y^TA >= 0
//    y^Tb <  0
// In the following instance, the first system has no solution,
// since adding up the two inequalities gives x_2 <= -1:
//     x_1 - 2x_2  <=  1
//    -x_1 + 3x_2  <= -2
//     x_1,   x_2  >=  0

int main() {
  int  Ax1[] = { 1, -1};                        // column for x1
  int  Ax2[] = {-2,  3};                        // column for x2
  int*   A[] = {Ax1, Ax2};                      // A comes columnwise
  int    b[] = {1, -2};                         // right-hand side
  CGAL::Comparison_result
    r[] = {CGAL::SMALLER, CGAL::SMALLER};      // constraints are "<="
  int    c[] = {0, 0};                         // zero objective function

  // now construct the linear program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program lp (2, 2, A, b, r, c);

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_linear_program(lp, ET());

  // get certificate for infeasibility
  assert (s.is_infeasible());
  Solution::Infeasibility_certificate_iterator y =
    s.infeasibility_certificate_begin();
  // check y >= 0
  assert (ET(y[0]) >= 0);
  assert (ET(y[1]) >= 0);
  // check y^T A >= 0
  assert (ET(y[0]) * A[0][0] + ET(y[1]) * A[0][1] >= 0);
  assert (ET(y[0]) * A[1][0] + ET(y[1]) * A[1][1] >= 0);
  // check y^T b < 0
  assert (ET(y[0]) * b[0] + ET(y[1]) * b[1] < 0);

  return 0;
}
