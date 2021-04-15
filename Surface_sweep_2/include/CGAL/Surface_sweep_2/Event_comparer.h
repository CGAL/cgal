// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_EVENT_COMPARER_H
#define CGAL_SURFACE_SWEEP_2_EVENT_COMPARER_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Comparison functor of curves used by the surface-sweep algorithm.
 */

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class
 * A functor used to compare events and event points in an xy-lexicographic
 * order. Used to maintain the order of the event queue (the X-structure)
 * in the sweep-line algorithm.
 */
template <typename GeometryTraits_2, typename Event_>
class Event_comparer {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Event_                                        Event;

private:
  typedef Geometry_traits_2                             Gt2;
public:

  typedef typename Gt2::Point_2                         Point_2;
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;

  // should be ok, as Gt2 is supposed to be the adaptor
  typedef typename Gt2::Left_side_category              Left_side_category;
  typedef typename Gt2::Bottom_side_category            Bottom_side_category;
  typedef typename Gt2::Top_side_category               Top_side_category;
  typedef typename Gt2::Right_side_category             Right_side_category;

private:
  // Data members:
  const Gt2* m_traits;                  // The geometric-traits object.

  Arr_parameter_space m_ps_in_x;        // Storing curve information when
  Arr_parameter_space m_ps_in_y;        // comparing a curve end with
  Arr_curve_end m_index;                // boundary conditions.

public:
  /*! Cosntruct. */
  Event_comparer(const Gt2* traits) : m_traits(traits) {}

  /*! Compare two existing events.
   * This operator is called by the multiset assertions only in
   * debug mode (to verify that event was inserted at the right place).
   */
  Comparison_result operator()(const Event* e1, const Event* e2) const
  {
    const bool on_boundary1 = e1->is_on_boundary();
    const bool on_boundary2 = e2->is_on_boundary();

    if (! on_boundary1 && ! on_boundary2) {
      // Both events do not have boundary conditions - just compare the points.
      return m_traits->compare_xy_2_object()(e1->point(), e2->point());
    }

    if (! on_boundary1) {
      // Compare the point associated with the first event with the second
      // boundary event.
      return this->operator()(e1->point(), e2);
    }

    if (! on_boundary2) {
      // Compare the point associated with the second event with the first
      // boundary event.
      return CGAL::opposite(this->operator()(e2->point(), e1));
    }

    return _compare_curve_end_with_event(e1->curve(), _curve_end(e1),
                                         e1->parameter_space_in_x(),
                                         e1->parameter_space_in_y(),
                                         e2);
  }

  /*! Compare a point, which should be inserted into the event queue,
   * with an existing event point.
   */
  Comparison_result operator()(const Point_2& pt, const Event* e2) const
  {
    const bool on_boundary2 = e2->is_on_boundary();

    if (! on_boundary2) {
      // If e2 is a normal event, just compare pt and the event point.
      return m_traits->compare_xy_2_object() (pt, e2->point());
    }

    // Get the sign of the event's boundary condition in x. Note that a valid
    // point is always larger than any negative boundary event and smaller
    // than any positive boundary event.
    Arr_parameter_space ps_x2 = e2->parameter_space_in_x();

    if (ps_x2 == ARR_LEFT_BOUNDARY) return LARGER;
    else if (ps_x2 == ARR_RIGHT_BOUNDARY) return SMALLER;

    // Get the curve end that e2 represents, and compare the x-position of the
    // given point and this curve end.
    Arr_curve_end ind = _curve_end(e2);
    Comparison_result res =
      m_traits->compare_x_point_curve_end_2_object()(pt, e2->curve(), ind);

    if (res != EQUAL) return res;

    // The event and the point has the same x-position. Get the sign of the
    // event's boundary condition in y. Note that a valid point is always
    // larger than any negative boundary event and smaller than any positive
    // boundary event.
    Arr_parameter_space ps_y2 = e2->parameter_space_in_y();

    CGAL_assertion (ps_y2 != ARR_INTERIOR);
    return (ps_y2 == ARR_BOTTOM_BOUNDARY) ? LARGER : SMALLER;
  }

  /*! Compare a curve end, which should be inserted into the event queue,
   * with an existing event point.
   * Note that the index of the curve end as well as its boundary conditions
   * must be set beforehand using set_index() and set_parameter_space_in_x/y().
   */
  Comparison_result operator()(const X_monotone_curve_2& cv,
                               const Event* e2) const
  {
    return _compare_curve_end_with_event(cv, m_index, m_ps_in_x, m_ps_in_y, e2);
  }

  /// \name Set the boundary conditions of a curve end we are about to compare.
  //@{
  void set_parameter_space_in_x(Arr_parameter_space bx) { m_ps_in_x = bx; }

  void set_parameter_space_in_y(Arr_parameter_space by) { m_ps_in_y = by; }

  void set_index(Arr_curve_end ind) { m_index = ind; }
  //@}

private:
  /*! Compare a given curve end with an event.
   * \param cv The curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \param e2 The event, which may have boundary conditions.
   * \return The comparison result of the curve end with the event.
   */
  Comparison_result
  _compare_curve_end_with_event(const X_monotone_curve_2& cv,
                                Arr_curve_end ind,
                                Arr_parameter_space ps_x,
                                Arr_parameter_space ps_y,
                                const Event* e2) const
  {
    // Check if the curve end has a boundary condition in x.
    if (ps_x == ARR_LEFT_BOUNDARY) {
      if (e2->parameter_space_in_x() == ARR_LEFT_BOUNDARY) {
        // Both defined on the left boundary - compare them there.
        CGAL_assertion (ind == ARR_MIN_END);

        return m_traits->compare_y_curve_ends_2_object() (cv, e2->curve(), ind);
      }

      // The curve end is obviously smaller.
      return SMALLER;
    }

    if (ps_x == ARR_RIGHT_BOUNDARY) {
      if (e2->parameter_space_in_x() == ARR_RIGHT_BOUNDARY) {
        // Both defined on the right boundary - compare them there.
        CGAL_assertion (ind == ARR_MAX_END);

        return m_traits->compare_y_curve_ends_2_object()(cv, e2->curve(), ind);
      }

      // The curve end is obviously larger.
      return LARGER;
    }

    // Check if the event has a boundary condition in x. Note that if it
    // has a negative boundary condition, the curve end is larger than it,
    // and if it has a positive boundary condition, the curve end is smaller.
    if (e2->parameter_space_in_x() == ARR_LEFT_BOUNDARY) return (LARGER);
    if (e2->parameter_space_in_x() == ARR_RIGHT_BOUNDARY) return (SMALLER);

    CGAL_assertion(ps_y != ARR_INTERIOR);
    Comparison_result res;

    Arr_curve_end ind2 = _curve_end(e2);

    // Act according to the boundary sign of the event.
    if (e2->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) {

      // Compare the x-positions of the two entities.
      res = m_traits->compare_x_curve_ends_2_object()(cv, ind,
                                                      e2->curve(), ind2);
      if (res != EQUAL) return res;

      // In case of equal x-positions, the curve end is larger than the event,
      // which lies on the bottom boundary (unless it also lies on the bottom
      // boundary).
      if (ps_y == ARR_BOTTOM_BOUNDARY) return EQUAL;

      return (LARGER);
    }

    if (e2->parameter_space_in_y() == ARR_TOP_BOUNDARY) {
      // Compare the x-positions of the two entities.
      res =
        m_traits->compare_x_curve_ends_2_object()(cv, ind, e2->curve(), ind2);
      if (res != EQUAL) return res;

      // In case of equal x-positions, the curve end is smaller than the event,
      // which lies on the top boundary (unless it also lies on the top
      // boundary).
      if (ps_y == ARR_TOP_BOUNDARY) return EQUAL;

      return SMALLER;
    }

    // If we reached here, e2 is not a boundary event and is associated with
    // a valid point. We compare the x-position of this point with the curve
    // end.
    res = m_traits->compare_x_point_curve_end_2_object()(e2->point(), cv, ind);

    if (res != EQUAL) return CGAL::opposite(res);

    // In case of equal x-positions, is the curve end has a negative boundary
    // sign, then it lies on the bottom boundary below the event. Otherwise,
    // it lies on the top aboundary above the event e2.
    return (ps_y == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
  }

  /*! Detemine if the given event represents a left or a right curve end. */
  inline Arr_curve_end _curve_end(const Event* e) const
  {
    return (e->has_left_curves()) ?
      ((e->is_right_end()) ? ARR_MAX_END : ARR_MIN_END) :
      ((e->is_left_end()) ? ARR_MIN_END : ARR_MAX_END);
  }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
