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

#include <CGAL/QP_solver/gmp_double.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_full_exact_pricing.h>
#include <CGAL/QP_partial_exact_pricing.h>
#include <CGAL/QP_full_filtered_pricing.h>
#include <CGAL/QP_partial_filtered_pricing.h>
#include <CGAL/QP_solver/Double.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>

#include <CGAL/QP_solver/MPS.h> // should to into QP_solver.h (?)

int main(const int argNr,const char **args) {
  using std::cout;
  using std::endl;
  using CGAL::Tag_true;
  using CGAL::Tag_false;

  // get desired level of additional logging output:
  const int verbosity = argNr < 2? 1 : std::atoi(args[1]);

  // construct QP instance:
  typedef double IT;
  typedef CGAL::Double ET;
  typedef CGAL::QP_MPS_instance<IT,ET> QP;
  QP qp(std::cin,true,verbosity);

  // check for format errors in MPS file:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
         << "Error: " << qp.error() << endl;
    std::exit(2);
  }

  if (verbosity > 0) {
    cout << endl << qp << endl;
  }

  typedef Tag_false Is_linear; // is the instance known in advance to be an LP?
  typedef Tag_false Is_symmetric; // is the D matrix known to be symmetric?
  typedef Tag_false Has_equalities_only_and_full_rank; // (see manual)
  typedef Tag_false Is_in_standard_form; // (see manual)?

  // in case of an LP, zero the D matrix:
  // (Note: if you know in advance that the problem is an LP
  // you should not do this, but set Is_linear to Tag_true.)
  if (qp.is_linear() && !check_tag(Is_linear()))
    qp.make_zero_D();

  typedef CGAL::QP_solver_MPS_traits_d<QP,
    Is_linear,
    Is_symmetric,
    Has_equalities_only_and_full_rank,
    Is_in_standard_form>         Traits;

  CGAL::QP_pricing_strategy<Traits> *strategy = 0;
  //    new CGAL::QP_partial_filtered_pricing<Traits,IT>;

  typedef CGAL::QP_solver<Traits> Solver;
  Solver solver(qp.number_of_variables(),
                                 qp.number_of_constraints(),
                                 qp.A(),qp.b(),qp.c(),
                                 qp.D(),
                                 qp.row_types(),
                                 qp.fl(),qp.l(),qp.fu(),qp.u(),
                                 strategy,
                                 verbosity);

  if (solver.is_valid()) {
    cout << "Solution is valid." << endl;
  } else {
    cout << "Solution is not valid!" << endl;
    return 1;
  }

  // get solution:
  if (solver.status() != Solver::INFEASIBLE) {
    cout << "Solution:" << endl;
    std::vector<ET> x(qp.number_of_variables(),0);  // solution vector
    // Note: following should be fixed to also output values for nonbasics...
    Solver::Basic_variable_numerator_iterator
      x_it =  solver.basic_original_variables_numerator_begin(),
      x_end = solver.basic_original_variables_numerator_end();
    Solver::Basic_variable_index_iterator 
      x_ind = solver.basic_original_variables_index_begin();
    const ET denom = solver.variables_common_denominator();
    for (; x_it != x_end; ++x_it, ++x_ind)
      x[*x_ind] = *x_it;
    
    // output solution:
    cout << "Exact solution:" << endl;
    for (unsigned int i=0; i<qp.number_of_variables(); ++i)
      cout << "x[" << i << "] = " << x[i] << "/" << denom << endl;
    cout << "Solution (as doubles):" << endl;
  for (unsigned int i=0; i<qp.number_of_variables(); ++i)
    cout << "x[" << i << "] ~= "
         << CGAL::to_double(x[i])/CGAL::to_double(denom)
         << endl;
  }
}
