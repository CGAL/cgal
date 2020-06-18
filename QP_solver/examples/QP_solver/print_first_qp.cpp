#include <iostream>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// program type
typedef CGAL::Quadratic_program<int> Program;

int main() {
  // by default, we have a nonnegative QP with Ax <= b
  Program qp (CGAL::SMALLER, true, 0, false, 0);

  // now set the non-default entries: 0 <-> x, 1 <-> y
  qp.set_a(0, 0,  1); qp.set_a(1, 0, 1); qp.set_b(0, 7);  //  x + y  <= 7
  qp.set_a(0, 1, -1); qp.set_a(1, 1, 2); qp.set_b(1, 4);  // -x + 2y <= 4
  qp.set_u(1, true, 4);                                   //       y <= 4
  qp.set_d(0, 0, 2); qp.set_d (1, 1, 8);                  // x^2 + 4 y^2
  qp.set_c(1, -32);                                       // -32y
  qp.set_c0(64);                                          // +64

  // print the program in MPS format
  CGAL::print_nonnegative_quadratic_program
    (std::cout, qp, "first_qp");

  return 0;
}
