// example: read linear program in MPS format from file
// the LP below is the first linear program example in the user manual
#include <iostream>
#include <fstream>

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
  std::cout << s;
  return 0;
}
