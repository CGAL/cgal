#include <iostream>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// program type
typedef CGAL::Quadratic_program<int> Program;

int main() {
  // by default, we have a nonnegative QP with Ax <= b
  Program lp (CGAL::SMALLER, true, 0, false, 0);

  // now set the non-default entries: 0 <-> x, 1 <-> y
  lp.set_a(0, 0,  1); lp.set_a(1, 0, 1); lp.set_b(0, 7);  //  x + y  <= 7
  lp.set_a(0, 1, -1); lp.set_a(1, 1, 2); lp.set_b(1, 4);  // -x + 2y <= 4
  lp.set_c(1, -32);                                       // -32y

  // print the program in MPS format
  CGAL::print_nonnegative_linear_program
    (std::cout, lp, "first_nonnegative_lp");

  return 0;
}
