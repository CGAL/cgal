// example: invert a random matrix (or find out that it's singular)
#include <iostream>
#include <vector>
#include <algorithm>

#include <CGAL/Random.h>
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

int n = 4;         // dimension of matrix
int s = 10;        // coordinates are in {-s, -s+1,...,0,...,s-1,s}
CGAL::Random rd;   // random number generator

// program and solution types
typedef CGAL::Quadratic_program<int> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

int main() {
  std::vector<std::vector<CGAL::Quotient<ET> > >
    inv_a; // stored columnwise

  // we need a free LP (no variable bounds) with Ax = b
  Program lp (CGAL::EQUAL, false, 0, false, 0);

  // constraint matrix A: the random matrix to be inverted
  for (int j=0; j<n; ++j)
    for (int i=0; i<n; ++i)
      lp.set_a (j, i, rd.get_int (-s, s));

  // we need to solve n LP, one for every right-hand side b =  e_j to
  // get j-th column of the inverse
  for (int j=0; j<n; ++j) {
    lp.set_b (j, 1);

    // solve the lp, using ET as the exact type
    Solution s = CGAL::solve_linear_program(lp, ET());
    if (s.is_infeasible()) {
      std::cout << "matrix is singular" << std::endl;
      return 0;
    } else {
      // store solution
      inv_a.push_back(std::vector<CGAL::Quotient<ET> >());
      std::copy (s.variable_values_begin(), s.variable_values_end(),
                 std::back_inserter (inv_a[j]));
    }
    lp.set_b (j, 0); // reset for next round
  }

  // output both matrices, and check that they are indeed inverse to each
  // other
  std::cout << "Random matrix A...:" << std::endl;
  Program::A_iterator a = lp.get_a();
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      std::cout << (*(a+j))[i] << " "; // row i
    std::cout << std::endl;
  }

  std::cout << std::endl << "...and its inverse: " << std::endl;
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      std::cout << inv_a[j][i] << " "; // row i
    std::cout << std::endl;
  }

  // check inverse property
  for (int i=0; i<n; ++i)
    for (int j=0; j<n; ++j) {
      // i-th row of A times j-th column of inverse
      CGAL::Quotient<ET> val;
      for (int k=0; k<n; ++k) val += (*(a+k))[i] * inv_a[j][k];
      assert (val == (i == j ? 1 : 0));
    }

  return 0;
}
