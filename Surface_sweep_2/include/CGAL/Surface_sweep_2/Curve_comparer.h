// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_CURVE_COMPARER_H
#define CGAL_SURFACE_SWEEP_2_CURVE_COMPARER_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Comparison functor of curves used by the surface-sweep algorithm.
 */

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! A functor used to compare curves and curve endpoints by their vertical
 * y-order. Used to maintain the order of curves in the status line (the
 * Y-structure) in the surface-sweep algorithm.
 */
template <typename GeometryTraits_2, typename Event_, typename Subcurve_>
class Curve_comparer {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef Arr_traits_basic_adaptor_2<Gt2>               Traits_adaptor_2;

  typedef typename Traits_adaptor_2::Point_2            Point_2;
  typedef typename Traits_adaptor_2::X_monotone_curve_2 X_monotone_curve_2;

private:
  const Traits_adaptor_2* m_traits;     // A geometric-traits object.
  Event** m_curr_event;                 // Points to the current event point.

public:
  /*! Construct. */
  Curve_comparer(const Traits_adaptor_2* t, Event** e_ptr) :
    m_traits(t),
    m_curr_event(e_ptr)
  {}

  /*! Compare the vertical position of two subcurves in the status line.
   * This operator is called only in debug mode.
   */
  Comparison_result operator()(const Subcurve* c1, const Subcurve* c2) const
  {
    // In case to two curves are right curves at the same event, compare
    // to the right of the event point.
    if (std::find((*m_curr_event)->right_curves_begin(),
                  (*m_curr_event)->right_curves_end(),
                  c1) != (*m_curr_event)->right_curves_end() &&
        std::find((*m_curr_event)->right_curves_begin(),
                  (*m_curr_event)->right_curves_end(),
                  c2) != (*m_curr_event)->right_curves_end())
    {
      return m_traits->compare_y_at_x_right_2_object()
              (c1->last_curve(), c2->last_curve(), (*m_curr_event)->point());
    }

    Arr_parameter_space ps_x1 =
      m_traits->parameter_space_in_x_2_object()(c1->last_curve(), ARR_MIN_END);
    Arr_parameter_space ps_y1 =
      m_traits->parameter_space_in_y_2_object()(c1->last_curve(), ARR_MIN_END);

    if ((ps_x1 == ARR_INTERIOR) && (ps_y1 == ARR_INTERIOR))
    {
      // The first curve has a valid left endpoint. Compare the y-position
      // of this endpoint to the second subcurve.
      return m_traits->compare_y_at_x_2_object()
        (m_traits->construct_min_vertex_2_object()(c1->last_curve()),
         c2->last_curve());
    }

    // We use the fact that the two curves are interior disjoint. As c2 is
    // already in the status line, then if c1 left end has a negative boundary
    // condition it obviously above it.
    CGAL_assertion(ps_x1 != ARR_RIGHT_BOUNDARY);

    if (ps_x1 == ARR_LEFT_BOUNDARY) return LARGER;

    // For similar reasons, if c1 begins on the bottom boundary it is below
    // c2, if it is on the top boundary it is above it.
    CGAL_assertion (ps_y1 != ARR_INTERIOR);
    return (ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
  }

  /*! Compare the relative y-order of the given point and the given subcurve.
   */
  Comparison_result operator()(const Point_2& pt, const Subcurve* sc) const
  { return m_traits->compare_y_at_x_2_object()(pt, sc->last_curve()); }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
