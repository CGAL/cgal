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

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <CGAL/basic.h>

#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
#else
#include <CGAL/Gmpz.h>
#endif

#include <CGAL/QP_solver.h>
#include <CGAL/QP_solver/QP_full_exact_pricing.h>
#include <CGAL/QP_solver/QP_exact_bland_pricing.h>
#include <CGAL/QP_solver/QP_partial_exact_pricing.h>
#include <CGAL/QP_solver/QP_full_filtered_pricing.h>
#include <CGAL/QP_solver/QP_partial_filtered_pricing.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_solution.h>

int main(const int argNr,const char **args) {
  using std::cout;
  using std::endl;
  using CGAL::Tag_true;
  using CGAL::Tag_false;

  // get desired level of additional logging output:
  const int verbosity = argNr < 2? 1 : std::atoi(args[1]);

  // construct QP instance:
  typedef int IT;
#ifndef CGAL_USE_GMP
  typedef CGAL::MP_Float ET;
#else
  typedef CGAL::Gmpz ET;
#endif
  typedef CGAL::QP_from_mps<IT, CGAL::Tag_false, CGAL::Tag_true> QP;
  QP qp(std::cin,true,verbosity);

  // check for format errors in MPS f\ile:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
         << "Error: " << qp.error() << endl;
    std::exit(2);
  }

  if (verbosity > 0) {
    CGAL::print_quadratic_program (cout, qp);
    cout << std::endl;
  }
  typedef CGAL::QP_solver_impl::QP_tags<CGAL::Tag_false,CGAL::Tag_false> Tags;
  CGAL::QP_pricing_strategy<QP, ET, Tags> *pricing_strategy =
    new CGAL::QP_full_exact_pricing<QP, ET, Tags>;
  typedef CGAL::QP_solver<QP, ET, Tags> Solver;
  Solver s (qp, pricing_strategy, 0);

  // get solution:
  if (s.status() == CGAL::QP_OPTIMAL) {
    // output solution:
    cout << "Objective function value: " << 
      s.solution() << endl;     
     
    cout << "Variable values:" << endl;
    Solver::Variable_value_iterator it 
      = s.original_variable_values_begin() ;
    for (int i=0; i < qp.n(); ++it, ++i)
      cout << "  " << qp.name_of_variable(i) << " = "
	   << CGAL::to_double(*it) << endl;
    return 0;
  }
  if  (s.status() == CGAL::QP_INFEASIBLE)
    cout << "Problem is infeasible." << endl;
  else // unbounded
    cout << "Problem is unbounded." << endl;

  delete pricing_strategy;
  return 0;
}
