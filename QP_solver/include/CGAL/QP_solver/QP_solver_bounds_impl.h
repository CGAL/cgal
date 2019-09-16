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

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::has_finite_lower_bound(int i) const
  // Given an index of an original or slack variable, returns whether
  // or not the variable has a finite lower bound.
{
  CGAL_qpe_assertion(i < qp_n + static_cast<int>(slack_A.size()));
  return i>=qp_n || check_tag(Is_nonnegative()) || *(qp_fl+i);
}

template < typename Q, typename ET, typename Tags >
bool QP_solver<Q, ET, Tags>::has_finite_upper_bound(int i) const
  // Given an index of an original or slack variable, returns whether
  // or not the variable has a finite upper bound.
{
  CGAL_qpe_assertion(i < qp_n + static_cast<int>(slack_A.size()));
  return i<qp_n && !check_tag(Is_nonnegative()) && *(qp_fu+i);
}

template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::lower_bound(int i) const
  // Given an index of an original or slack variable, returns its
  // lower bound.
{
  CGAL_qpe_assertion(i < qp_n + static_cast<int>(slack_A.size()));
  if (i < qp_n)                     // original variable?
    if (check_tag(Is_nonnegative()))
      return et0;
    else {
      CGAL_qpe_assertion(has_finite_lower_bound(i));
      return *(qp_l+i);
    }
  else                              // slack variable?
    return et0;
}

template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::upper_bound(int i) const
  // Given an index of an original variable, returns its upper bound.
{
  CGAL_qpe_assertion(i < qp_n); // Note: slack variables cannot have
				// finite upper bounds.
  CGAL_qpe_assertion(has_finite_upper_bound(i));
  return *(qp_u+i);
}

template < typename Q, typename ET, typename Tags >
typename QP_solver<Q, ET, Tags>::Bnd
QP_solver<Q, ET, Tags>::lower_bnd(int i) const
  // Given an index of an original, slack, or artificial variable,
  // return its lower bound.
{
  if (i < qp_n) {                                      // original variable?
    const bool is_finite = has_finite_lower_bound(i);
    return Bnd(false, is_finite, is_finite? lower_bound(i) : ET(0));
  } else                                              // slacky or art. var.?
    return Bnd(false, true, ET(0));
}

template < typename Q, typename ET, typename Tags >
typename QP_solver<Q, ET, Tags>::Bnd
QP_solver<Q, ET, Tags>::upper_bnd(int i) const
  // Given an index of an original, slack, or artificial variable,
  // return its upper bound.
{
  if (i < qp_n) {                                      // original variable?
    const bool is_finite = has_finite_upper_bound(i);
    return Bnd(true, is_finite, is_finite? upper_bound(i) : ET(0));
  } else                                              // slacky or art. var.?
    return Bnd(true, false, ET(0));
}

} //namespace CGAL

// ===== EOF ==================================================================
