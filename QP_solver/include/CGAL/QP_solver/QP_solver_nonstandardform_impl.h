// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer

namespace CGAL {

// Looks in x_O_v_i which bound is present for variable i and returns
// the variable's value corresponding to this bound.
//
// Precondition: Is_nonnegative is Tag_false.
template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::original_variable_value_under_bounds(int i) const
{
  CGAL_assertion(!check_tag(Is_nonnegative()) && i<qp_n);
  switch (x_O_v_i[i]) {
  case UPPER:
    return *(qp_u+i);
  case ZERO:
    return et0;
  case LOWER:
  case FIXED:
    return *(qp_l+i);
  case BASIC:
    CGAL_qpe_assertion(false);
  }
  return et0; // dummy
}

template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::variable_numerator_value(int i) const
{
  // Returns the current value of an *original* variable.
  CGAL_qpe_assertion( 0 <= i && i < qp_n );
  if (check_tag(Is_nonnegative())) {
    if (in_B[i] < 0) 
      return et0;
    else 
      return x_B_O[in_B[i]];
  }

  // now we have nonstandard form
  typedef QP_solver<Q, ET, Tags> QP;
  switch (x_O_v_i[i]) {
  case QP::UPPER:
    return ET(*(qp_u+i)) * d;
  case QP::ZERO:
    return et0;
  case QP::LOWER:
  case QP::FIXED:
    return ET(*(qp_l+i)) * d;
  case QP::BASIC:
    return x_B_O[in_B[i]];
  default: // never reached
    return et0;
  }
}

template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::nonbasic_original_variable_value
(int i) const
{
  if (check_tag(Is_nonnegative()))
    return et0;

  CGAL_assertion(!is_basic(i));
  return original_variable_value_under_bounds(i);
}

// Computes r_i:= A_i x_init, for i=row, where x_init is the solution
// with which the solver starts the computation. I.e., computes the
// scalar product of the row-th row of A and the vector x_init which
// contains as its entries the values original_variable_value(i),
// 0<=i<qp_n.
template < typename Q, typename ET, typename Tags >
ET  QP_solver<Q, ET, Tags>::multiply__A_ixO(int row) const
{
  ET value = et0;

  for (int i = 0; i < qp_n; ++i)
    // Note: the following computes
    //
    //   value += original_variable_value(i) * qp_A[i][row];
    //
    // but for efficiency, we only add summands that are known to be
    // nonzero.
    switch (x_O_v_i[i]) {
    case UPPER:
      value += ET(*(qp_u+i)) * ET(*((*(qp_A+i))+row));
      break;
    case LOWER:
    case FIXED:
      value += ET(*(qp_l+i)) * ET(*((*(qp_A+i))+row));
      break;
    case BASIC:
      CGAL_qpe_assertion(false);
    default:
      break;
    }

  return value;
}

// Computes r_{C}:= A_{C, N_O} x_{N_O}.
//
// Precondition: this routine should only be called for nonstandard form
// problems.
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__A_CxN_O(Value_iterator out) const
{
  CGAL_qpe_assertion(!check_tag(Is_nonnegative()));
  
  // initialize with zero vector:
  std::fill_n(out, C.size(), et0);
  
  for (int i = 0; i < qp_n; ++i)
    if (!is_basic(i)) {
      const ET x_i = nonbasic_original_variable_value(i);
      const A_column a_col = *(qp_A+i);
      Value_iterator out_it = out;
      for (Index_const_iterator row_it = C.begin();
	   row_it != C.end();
	   ++row_it, ++out_it)
	*out_it += x_i * ET(*(a_col+ *row_it));
    }
}

// Computes w:= 2D_{O, N_O} x_{N_O}.
//
// Precondition: this routine should only be called for nonstandard form
// problems.
//
// todo: In order to optimize this routine, we can if D is symmetric,
// multiply by two at the end of the computation instead of at each
// access to D. (Maybe its also faster to call
// nonbasic_original_variable_value() only O(n) times and not O(n^2)
// times.)
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__2D_OxN_O(Value_iterator out) const
{
  CGAL_qpe_assertion(!check_tag(Is_nonnegative()));

  // initialize with zero vector:
  std::fill_n(out, B_O.size(), et0);
  
  for (int row_it = 0; row_it < qp_n; ++row_it, ++out) {
    D_pairwise_accessor d_row(qp_D, row_it);
    for (int i = 0; i < qp_n; ++i)
      if (!is_basic(i)) {
	const ET value = nonbasic_original_variable_value(i);
	*out += d_row(i) * value;
      }
  }
}

// Computes r_{S_B}:= A_{S_B, N_O} x_{N_O}.
//
// Precondition: this routine should only be called for nonstandard form
// problems.
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__A_S_BxN_O(Value_iterator out) const
{
  // initialize with zero vector:
  std::fill_n(out, S_B.size(), et0);
  
  for (int i = 0; i < qp_n; ++i)
    if (!is_basic(i)) {
      const ET x_i = nonbasic_original_variable_value(i);
      const A_column a_col = *(qp_A+i);
      Value_iterator out_it = out;
      for (Index_const_iterator row_it = S_B.begin();
	   row_it != S_B.end();
	   ++row_it, ++out_it)
	*out_it += x_i * ET(*(a_col+ *row_it));
    }
}

// Initialize r_B_O.
//
// Note: this routine is called from transition() (and not during the
// initialization of the QP-solver).
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
init_r_B_O()
{
  CGAL_qpe_assertion(!check_tag(Is_nonnegative()) &&
			!check_tag(Is_linear()));
  r_B_O.resize(B_O.size());
  multiply__2D_B_OxN_O(r_B_O.begin());
}

// Initialize w.
//
// Note: this routine is called from transition() (and not during the
// initialization of the QP-solver).
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
init_w()
{
  CGAL_qpe_assertion(!check_tag(Is_nonnegative()) &&
			!check_tag(Is_linear()));
  w.resize(qp_n);
  multiply__2D_OxN_O(w.begin());
}

} //namespace CGAL

// ===== EOF ==================================================================
