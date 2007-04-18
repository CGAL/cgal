// example: read quadratic program in MPS format from file 
// the QP below is the first quadratic program example in the user manual
#include <iostream>
#include <fstream>
#include <CGAL/basic.h>
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
  std::ifstream in ("first_qp.mps");
  Program qp(in);         // read program from file
  assert (qp.is_valid()); // we should have a valid mps file

  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_quadratic_program(qp, ET());

  // output solution
  std::cout << s; 
  return 0;
}
