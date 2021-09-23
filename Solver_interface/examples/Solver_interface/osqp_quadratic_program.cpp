/*
* This example shows how to formulate and solve the following QP problem:
* https://osqp.org/docs/examples/setup-and-solve.html
*
*   Minimize
*   Objective:
*     1/2(4x1^2 + 2x1x2 + 2x2^2) + (x1 + x2) + 0 or in the matrix form
*     1/2 x^T|4 1|x + |1|^Tx + 0
*            |1 2|    |1|
*   Subject to
*     1 <= x1 + x2 <= 1
*     0 <= x1      <= 0.7
*     0 <=      x2 <= 0.7 or in the matrix form
*     |1|    |1 1|     |1.0|
*     |0| <= |1 0|x <= |0.7|
*     |0|    |0 1|     |0.7|
*
*   Expected results: x1 = 0.3; x2 = 0.7;
*/

#include <vector>
#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/OSQP_quadratic_program_traits.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;

int main(void) {

  const std::size_t n = 2; // number of variables
  const std::size_t m = 3; // number of constraints
  CGAL::OSQP_quadratic_program_traits<FT> osqp(n, m);

  osqp.set_P(0, 0, 4);
  osqp.set_P(0, 1, 1);
  osqp.set_P(1, 1, 2);

  osqp.set_q(0, 1);
  osqp.set_q(1, 1);

  osqp.set_r(0);

  osqp.set_A(0, 0, 1);
  osqp.set_A(0, 1, 1);
  osqp.set_A(1, 0, 1);
  osqp.set_A(2, 1, 1);

  osqp.set_l(0, 1);
  osqp.set_l(1, 0);
  osqp.set_l(2, 0);

  osqp.set_u(0, 1);
  osqp.set_u(1, 0.7);
  osqp.set_u(2, 0.7);

  std::vector<FT> x; x.reserve(2);
  osqp.solve(std::back_inserter(x));

  std::cout << "solution (x1 x2): ";
  for (const FT value : x) {
    std::cout << value << "; ";
  }
  std::cout << std::endl;
}
