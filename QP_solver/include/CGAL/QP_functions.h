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
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/QP_solver/include/CGAL/QP_functions.h $
// $Id: QP_functions.h 33922 2006-09-05 12:32:25Z gaertner $
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_QP_FUNCTIONS_H
#define CGAL_QP_FUNCTIONS_H

#include <CGAL/iterator.h>
#include <CGAL/QP_solver.h>
#include <CGAL/QP_models.h>

CGAL_BEGIN_NAMESPACE

template <class Q>
bool QP_is_in_equational_form (const Q& qp) {
  // check whether all constraints are equality constraints
  typename Q::R_iterator r = qp.r();
  typename Q::R_iterator r_end = r + qp.m();
  for (; r < r_end; ++r)
    if (*r != CGAL::EQUAL) return false;
  return true;
}

template <class Q, class ET>
bool QP_has_full_row_rank (const Q& qp, const ET& dummy) {
  // we solve the following auxiliary LP, using exact type ET:
  // --------
  // min 0
  // A x == 0
  //   x >= 0
  // --------
  // Then A has full row rank if and only if all artificials
  // have left the basis after phase I; the QP_solver diagnostics
  // tells us this
  //
  // auxiliary LP type
  typedef Const_oneset_iterator<typename Q::value_type>  C_iterator;
  typedef Const_oneset_iterator<typename Q::value_type>  B_iterator;
  typedef Const_oneset_iterator<CGAL::Comparison_result> R_iterator;
  typedef Nonnegative_LP_from_iterators
    <typename Q::A_iterator, B_iterator, R_iterator, C_iterator> LP;

  //  auxiliary LP
  LP lp (qp.n(), qp.m(), 
	 qp.a(), B_iterator(0), R_iterator(CGAL::EQUAL), C_iterator(0));

  //  solver Tags
  typedef QP_solver_impl::QP_tags<
    Tag_true,  // Is_linear
    Tag_true,  // Is_symmetric
    Tag_false, // Has_equalities_only_and_full_rank (we don't know it yet)
    Tag_true>  // Is_in_standard_form
  Tags;

  // solver type
  typedef QP_solver<LP, ET, Tags> Solver;

  // now solve auxiliary LP and compute predicate value
  Solver solver (lp);
  return !solver.diagnostics.redundant_equations;
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_FUNCTIONS_H
