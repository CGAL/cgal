// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
#include <CGAL/basic.h>
#include <iostream>
#include <sstream>
#include <cstdlib>

#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
#else
#include <CGAL/Gmpzf.h>
#endif 

#include <CGAL/QP_solver.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>


int main(const int argNr,const char **args) {
  using std::cout;
  using std::endl;

  // get desired level of additional logging output:
  const int verbosity = argNr < 2? 1 : std::atoi(args[1]);

  // construct QP instance:
  typedef double IT;

#ifndef CGAL_USE_GMP 
  typedef CGAL::MP_Float ET; 
#else
  typedef CGAL::Gmpzf ET;
#endif

  typedef CGAL::QP_from_mps<IT, CGAL::Tag_false, CGAL::Tag_true, CGAL::Tag_true> QP;
  QP qp(std::cin,true,verbosity);
  // check for format errors in MPS f\ile:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
         << "Error: " << qp.error() << endl;
    std::exit(2);
  }

//   if (verbosity > 0) {
//     CGAL::print_quadratic_program (cout, qp);
//     cout << std::endl;
//   }

  typedef CGAL::QP_solution<ET> Solution;
  Solution s = CGAL::solve_quadratic_program (qp, ET(0));

  if (s.is_valid()) {
    cout << "Solution is valid." << endl;
  } else {
    cout << "Solution is not valid!" << endl;
    return 1;
  }

  // get solution:
  if (s.status() == CGAL::QP_OPTIMAL) {
    // output solution:
    cout << "Objective function value: " << 
      CGAL::to_double(s.solution()) << endl;     
     
    cout << "Variable values:" << endl;
    Solution::Variable_value_iterator vit = 
      s.variable_values_begin() ;
    for (unsigned int i=0; i < qp.n(); ++vit, ++i)
      cout << "  " << qp.name_of_variable(i) << " = "
	   << CGAL::to_double(*vit) << endl;

    cout << "Basic variables (index, value):" << endl;
    Solution::Index_iterator iit = 
      s.basic_variable_indices_begin();
    vit = s.variable_values_begin();
    for ( ; iit < s.basic_variable_indices_end(); ++iit)
      cout << "  (" << *iit << ", " 
	   <<  CGAL::to_double(vit[*iit]) << ")" << endl;
    cout << "  There are " << s.number_of_basic_variables() 
	 << " basic variables. " << endl;

    cout << "Basic constraints: " << endl;
    Solution::Index_iterator cit = 
      s.basic_constraint_indices_begin();
    for ( ; cit < s.basic_constraint_indices_end(); ++cit)
      cout << *cit << " ";
    cout << endl;
    cout << "  There are " << s.number_of_basic_constraints() 
	 << " basic constraints. " << endl;

    return 0;
  }
  if  (s.status() == CGAL::QP_INFEASIBLE)
    cout << "Problem is infeasible." << endl;
  else // unbounded
    cout << "Problem is unbounded." << endl;
  return 0;
 
}
