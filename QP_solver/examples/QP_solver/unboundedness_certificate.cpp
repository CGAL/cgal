// example: extracting and verifying a proof of unboundedness from the solution
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

// we demonstrate the unboundedness certificate: either the feasible
// linear program
// min c^T x
//       A x <= b
//         x >= 0
// is bounded, or there exists w such that
//         w >= 0
//       A w <= 0
//     c^t w <  0,
// in which case the objective function becomes arbitrarily small
// along the ray {x* + tw | t >= 0}, x* any feasible solution.
//
// In the following instance, the linear program is unbounded, since
// the ray {(t,t)| t>= 1} is feasible, and the objective function becomes
// arbitrarily small on it:
// min -x_1 - 2x_2
//      x_1 - 2x_2  <=  -1
//     -x_1 +  x_2  <=   2
//      x_1,   x_2  >=   0

int main() {
  int  Ax1[] = { 1, -1};                        // column for x1
  int  Ax2[] = {-2,  1};                        // column for x2
  int*   A[] = {Ax1, Ax2};                      // A comes columnwise
  int    b[] = {-1, 2};                         // right-hand side
  CGAL::Comparison_result
    r[] = {CGAL::SMALLER, CGAL::SMALLER};       // constraints are "<="
  int    c[] = {-1, -2};                        // objective function

  // now construct the linear program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program lp (2, 2, A, b, r, c);

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_nonnegative_linear_program(lp, ET());

  // get certificate for unboundedness
  assert (s.is_unbounded());
  Solution::Unboundedness_certificate_iterator w =
    s.unboundedness_certificate_begin();
  // check w >= 0
  assert (ET(w[0]) >= 0);
  assert (ET(w[1]) >= 0);
  // check A w <= 0
  assert (A[0][0] * ET(w[0]) + A[1][0] * ET(w[1]) <= 0);
  assert (A[0][1] * ET(w[0]) + A[1][1] * ET(w[1]) <= 0);
  // check c^T w < 0
  assert (c[0] * ET(w[0]) + c[1] * ET(w[1]) < 0);

  return 0;
}
