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
// implementation: test program for the QP solver
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

template <typename T>
T string_to(const std::string& s) {
  std::stringstream strm(s);
  T t;
  strm >> t;
  return t;
}

void bailout(const char *msg)
{
  std::cout << "Error: " << msg << '.' << std::endl;
  exit(1);
}

int main(const int argNr,const char **args) {
  using std::cout;
  using std::endl;
  using CGAL::Tag_true;
  using CGAL::Tag_false;

  // get desired level of additional logging output:
  const int verbosity = argNr < 2? 1 : string_to<int>(args[1]);

  // construct QP instance:
  typedef double IT;
  typedef CGAL::Double ET;
  typedef CGAL::QP_solver_MPS_traits_d<IT,ET,
    Tag_false, // is the instance known in advance to be an LP?
    Tag_false, // is the instance's D matrix known to be symmetric?
    Tag_false, // Has_equalities_only_and_full_rank (see manual)?
    Tag_true>  // Is_in_standard_form (see manual)?
    Traits;
  CGAL::QP_MPS_instance<Traits> qp(std::cin);

  // check for format errors in MPS file:
  if (!qp.is_valid()) {
    cout << "Input is not a valid MPS file." << endl
	 << "Error: " << qp.error() << endl;
    exit(2);
  }

  if (verbosity > 0) {
    cout << endl << qp << endl;
  }

  CGAL::QP_solver<Traits> solver(qp.number_of_variables(),
				 qp.number_of_constraints(),
				 qp.A(),qp.b(),qp.c(),qp.D(),
				 qp.row_types(),
				 qp.default_pricing_strategy());

  return 0;
}
