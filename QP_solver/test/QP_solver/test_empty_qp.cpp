#include <cstdlib>
#include <cassert>
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

int main() 
{
  CGAL::Quadratic_program<int> qp; // no variables, no constraints
  assert (qp.is_valid());

  CGAL::Quadratic_program_solution<ET> s = 
    CGAL::solve_quadratic_program (qp, ET());
  
  std::cout << "Solution = " << s.solution() << std::endl;
  std::cout << "Status   = " << s.status() << std::endl;

  return 0;
}
