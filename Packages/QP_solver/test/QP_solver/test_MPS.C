// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : test/QP_solver/test_MPS.C
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.2
// revision_date : 2000/08/21
//
// author(s)     : Kaspar Fischer (fischerk@inf.ethz.ch)
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: example QP solver for QP's in MPS format
// ============================================================================

#include <iostream>
#include <sstream>
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
    exit(2);
  }

  if (verbosity > 0) {
    cout << endl << qp << endl;
  }

  typedef Tag_false Is_linear; // is the instance known in advance to be an LP?
  typedef Tag_false
 Is_symmetric; // is the D matrix known to be symmetric?
  typedef Tag_false Has_equalities_only_and_full_rank; // (see manual)
  typedef Tag_true Is_in_standard_form; // (see manual)?

  // in case of an LP, zero the D matrix:
  // (Note: if you known in advance that the problem is an LP
  // you should not do this, see below.)
  if (qp.is_linear() && !check_tag(Is_linear()))
    qp.make_zero_D();

  typedef CGAL::QP_solver_MPS_traits_d<QP,
    Is_linear,
    Is_symmetric,
    Has_equalities_only_and_full_rank,
    Is_in_standard_form>         Traits;

  CGAL::QP_pricing_strategy<Traits> *strategy = 0;
  //    new CGAL::QP_full_filtered_pricing<Traits,IT>;

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
    cout << "Solution is not valid!." << endl;
    return 1;
  }

  // get solution:
  cout << "Solution:" << endl;
  std::vector<ET> x(qp.number_of_variables(),0);  // solution vector
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
