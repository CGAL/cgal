// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <CGAL/config.h>

#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
#else
#include <CGAL/Gmpz.h>
#endif

#include <CGAL/QP_solver/QP_solver.h>
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
  typedef CGAL::Quadratic_program_from_mps<IT> QP;
  QP qp(std::cin);

  // check for format errors in MPS f\ile:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
         << "Error: " << qp.get_error() << endl;
    std::exit(2);
  }

  if (verbosity > 0) {
    CGAL::print_quadratic_program (cout, qp);
    cout << std::endl;
  }

  CGAL::Quadratic_program_options options;
  options.set_pricing_strategy(CGAL::QP_DANTZIG);
  options.set_verbosity(0);
  typedef CGAL::Quadratic_program_solution<ET> Solution;
  Solution s = CGAL::solve_quadratic_program (qp, ET(), options);

  // get solution:
  if (s.status() == CGAL::QP_OPTIMAL) {
    // output solution:
    cout << "Objective function value: " <<
      s.objective_value() << endl;

    cout << "Variable values:" << endl;
    Solution::Variable_value_iterator it
      = s.variable_values_begin() ;
    for (int i=0; i < qp.get_n(); it++, ++i)
      cout << "  " << qp.variable_name_by_index(i) << " = "
           << CGAL::to_double(*it) << endl;
    return 0;
  }
  if  (s.status() == CGAL::QP_INFEASIBLE)
    cout << "Problem is infeasible." << endl;
  else // unbounded
    cout << "Problem is unbounded." << endl;
  return 0;
}
