// Copyright (c) 2006,2007,2009,2010,2011,2014 Tel-Aviv University (Israel).
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
//
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_SWEEP_LINE_FUNCTORS_H
#define CGAL_SWEEP_LINE_FUNCTORS_H

/*! \file
 * Comparison functors used by the sweep-line algorithm.
 */

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL {

/*! \class
 * A functor used to compare events and event points in an xy-lexicographic
 * order. Used to maintain the order of the event queue (the X-structure)
 * in the sweep-line algorithm.
 */
template <class Traits_, class Event_>
class Compare_events
{
public:

  typedef Traits_                                   Traits_2;
  typedef Event_                                    Event;

  typedef typename Traits_2::Point_2                Point_2;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;

  // should be ok, as Traits_2 is supposed to be the adaptor
  typedef typename Traits_2::Left_side_category          Left_side_category;
  typedef typename Traits_2::Bottom_side_category        Bottom_side_category;
  typedef typename Traits_2::Top_side_category           Top_side_category;
  typedef typename Traits_2::Right_side_category         Right_side_category;

private:

  // Data members:
  const Traits_2 * m_traits;                // The geometric-traits object.

  Arr_parameter_space  m_ps_in_x;           // Storing curve information when
  Arr_parameter_space  m_ps_in_y;           // comparing a curve end with
  Arr_curve_end        m_ind;               // boundary conditions.

public:

  /*! Cosntructor. */
  Compare_events (const Traits_2 * traits) :
    m_traits (traits)
  {}

  /*!
   * Compare two existing events.
   * This operator is called by the multiset assertions only in
   * debug mode (to verify that event was inserted at the right place).
   */
  Comparison_result operator()(const Event* e1, const Event* e2) const
  {
    const bool  on_boundary1 = e1->is_on_boundary();
    const bool  on_boundary2 = e2->is_on_boundary();

    if (! on_boundary1 && ! on_boundary2)
    {
      // Both events do not have boundary conditions - just compare the points.
      return (m_traits->compare_xy_2_object()(e1->point(), e2->point()));
    }

    if (! on_boundary1)
    {
      // Compare the point associated with the first event with the second
      // boundary event.
      return ( this->operator()(e1->point(), e2) );
    }

    if (! on_boundary2)
    {
      // Compare the point associated with the second event with the first
      // boundary event.
      return (CGAL::opposite(this->operator()(e2->point(), e1)));
    }

    return (_compare_curve_end_with_event (e1->curve(), _curve_end(e1),
                                           e1->parameter_space_in_x(),
                                           e1->parameter_space_in_y(),
                                           e2));
  }

  /*!
   * Compare a point, which should be inserted into the event queue,
   * with an existing event (which might lie on the boundary)
   */
  Comparison_result operator() (const Point_2& pt, const Event* e2) const
  {
    Arr_parameter_space ps_x1 = m_traits->parameter_space_in_x_2_object()(pt);
    Arr_parameter_space ps_y1 = m_traits->parameter_space_in_y_2_object()(pt);

    return _compare_point_with_event(pt, ps_x1, ps_y1, e2);
  }

  /*!
   * Compare a curve end, which should be inserted into the event queue,
   * with an existing event point.
   * Note that the ind of the curve end as well as its boundary conditions
   * must be set beforehand using set_ind() and set_parameter_space_in_x/y().
   */
  Comparison_result operator() (const X_monotone_curve_2& cv,
                                const Event* e2) const
  {
    return _compare_curve_end_with_event(cv, m_ind, m_ps_in_x, m_ps_in_y, e2);
  }

  /// \name Set the boundary conditions of a curve end we are about to compare.
  //@{
  void set_parameter_space_in_x (Arr_parameter_space bx)
  {
    m_ps_in_x = bx;
  }

  void set_parameter_space_in_y (Arr_parameter_space by)
  {
    m_ps_in_y = by;
  }

  void set_ind (Arr_curve_end ind)
  {
    m_ind = ind;
  }
  //@}

private:

  /*!
   * Compare a given point with an event. Note that \param pt
   * is the a curve-end ps_x1 and ps_y1 can have different values than
   * if \param pt would be an isolated point. The connected curve
   * determines ps_x1 and ps_y1 for a point that is a curve-end.
   *
   * \param pt The pint
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \param e2 The event, which may have boundary conditions.
   * \return The comparison result of the curve end with the event.
   */
  Comparison_result
  _compare_point_with_event (const Point_2& pt,
                             Arr_parameter_space ps_x1,
                             Arr_parameter_space ps_y1,
                             const Event* e2) const {

    std::cout << std::endl << "FUNCTOR pt: " << pt << std::endl;
    Arr_parameter_space ps_x2 = e2->parameter_space_in_x();
    Arr_parameter_space ps_y2 = e2->parameter_space_in_y();

    // e1 left      e2 left      compare_y_on_boundary(pt1, pt2) REMARK: obviously points on the left boundary are allowed
    // e1 left      e2 bottom    SMALLER
    // e1 left      e2 interior  SMALLER
    // e1 left      e2 top       SMALLER
    // e1 left      e2 right     SMALLER

    // e1 bottom    e2 left      LARGER
    // e1 bottom    e2 bottom    compare_x_on_boundary(pt1, pt2) REMARK: obviously points on the bottom boundary are allowed
    // e1 bottom    e2 interior  compare_x_on_boundary(pt1, pt2) REMARK: as one is on the boundary
    // e1 bottom    e2 top       compare_x_on_boundary(pt1, pt2) REMARK: as one is on the boundary
    // e1 bottom    e2 right     SMALLER

    // e1 interior  e2 left      LARGER
    // e1 interior  e2 bottom    compare_x_on_boundary(pt1, pt2) REMARK: as one is on the boundary
    // e1 interior  e2 interior  compare_xy(pt, pt)
    // e1 interior  e2 top       compare_x_on_boundary(pt1, pt2) REMARK: as one is on the boundary
    // e1 interior  e2 right     SMALLER

    // e1 top       e2 left      LARGER
    // e1 top       e2 bottom    compare_x_on_boundary(pt1, pt2) REMARK: as one is on the boundary
    // e1 top       e2 interior  compare_x_on_boundary(pt1, pt2) REMARK: as one is on the boundary
    // e1 top       e2 top       compare_x_on_boundary(pt1, pt2) REMARK: obviously points on the top boundary are allowed
    // e1 top       e2 right     SMALLER

    // e1 right     e2 left      LARGER
    // e1 right     e2 bottom    LARGER
    // e1 right     e2 interior  LARGER
    // e1 right     e2 top       LARGER
    // e1 right     e2 right     compare_y_on_boundary(pt1, pt2) REMARK: obviously points on the right boundary are allowed

    if (ps_x1 == ps_x2) {
      // same x-partition
      if (ps_x1 != ARR_INTERIOR) {
        // 1L2L, 1R2R
        // e2->point() is accessible, as pt is a point on the SAME left/right side
        Comparison_result res = (m_traits->compare_y_on_boundary_2_object() (pt, e2->point()));
        std::cout << "res1 " << res << std::endl;
        return res;
      }
      // else both are x-interior
      if ((ps_y1 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
        // both are y-interior, too 1I2I:
        // e2->point() is accessible as e2 is INTERIOR
        // NOTE compare_y_on_boundary could take a point that actually lies on the bottom or top side (if curve-end of a a curve on the identifaction)!!!!
        Comparison_result res = (m_traits->compare_xy_2_object() (pt, e2->point()));
        std::cout << "res2 " << res << std::endl;
        return res;
      } else {
        // at least onex2 of pt or e2 lies on a boundary
        // call the comparison depending on whether closed side is available
        typedef CGAL::internal::Arr_bottom_top_implementation_dispatch<
          Bottom_side_category, Top_side_category > BT;
        if (e2->is_closed()) {
          Comparison_result res = (_compare_x_point_with_closed_event(pt, e2, typename BT::Compare_x_on_boundary_2_points_tag()));
          std::cout << "res3 " << res << std::endl;
          return res;
        } else {
          Comparison_result res = (_compare_x_point_with_open_event  (pt, e2, typename BT::Compare_x_on_boundary_2_point_curve_end_tag()));
          std::cout << "res4 " << res << std::endl;
        return res;
        }
      }
    }

    // ps_x1 != ps_x2
    // only the 14 simple ps_x different cases remain:
    // 1L2B, 1L2I, 1L2T, 1L2R,
    // 1B2L, 1B2R,
    // 1I2L, 1I2R,
    // 1T2L, 1T2R
    // 1R2L, 1R2B, 1R2I, 1R2T
    if (ps_x1 == ARR_LEFT_BOUNDARY) {
      std::cout << "res5 SMALLER" << std::endl;
      return (SMALLER);
    }
    if (ps_x1 == ARR_RIGHT_BOUNDARY) {
      std::cout << "res6 LARGER" << std::endl;
      return (LARGER);
    }
    // ps_x1 == ARR_INTERIOR != ps_x2
    if (ps_x2 == ARR_LEFT_BOUNDARY) {
      std::cout << "res7 LARGER" << std::endl;
      return (LARGER);
    }
    // ps_x2 == ARR_RIGHT_BOUNDARY
    std::cout << "res8 SMALLER" << std::endl;
    return (SMALLER);
  }

  Comparison_result _compare_x_point_with_closed_event(const Point_2& pt, const Event* e2,
                                                     Arr_use_traits_tag) const {
    CGAL_assertion(e2->is_closed());
    std::cout << "cmpA"<< std::endl;
    return (m_traits->compare_x_on_boundary_2_object() (pt, e2->point()));
  }

  Comparison_result _compare_x_point_with_closed_event(const Point_2& pt, const Event* e2,
                                                       Arr_use_dummy_tag) const {
    CGAL_error();
    return (EQUAL);
  }

  Comparison_result _compare_x_point_with_open_event(const Point_2& pt, const Event* e2,
                                                     Arr_use_traits_tag) const {
    CGAL_assertion(!(e2->is_closed()));
    // e2 must be an open event here, and thus exactly one curve is incident
    std::cout << "cmpB"<< std::endl;
    return (m_traits->compare_x_point_curve_end_2_object() (pt, e2->curve(), _curve_end(e2)));
  }

  Comparison_result _compare_x_point_with_open_event(const Point_2& pt, const Event* e2,
                                                     Arr_use_dummy_tag) const {
    CGAL_error();
    return (EQUAL);
  }

  /*!
   * Compare a given curve end with an event.
   * \param cv The curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   * \param e2 The event, which may have boundary conditions.
   * \return The comparison result of the curve end with the event.
   */
  Comparison_result
  _compare_curve_end_with_event (const X_monotone_curve_2& cv,
                                 Arr_curve_end ind,
                                 Arr_parameter_space ps_x1,
                                 Arr_parameter_space ps_y1,
                                 const Event* e2) const
  {

    // if the curve-end is closed, use the point signature, which simplifies life
    if (m_traits->is_closed_2_object() (cv, ind)) {

      const Point_2& pt = (ind == ARR_MIN_END ?
                           m_traits->construct_min_vertex_2_object() (cv) :
                           m_traits->construct_max_vertex_2_object() (cv));

      // call point-comparsion, BUT use ps_x1, ps_y1 of curve-end!!!!
      return _compare_point_with_event(pt, ps_x1, ps_y1, e2);
    }

    // else cv has an OPEN end and ind, continue the "open way":

    Arr_parameter_space ps_x2 = e2->parameter_space_in_x();
    Arr_parameter_space ps_y2 = e2->parameter_space_in_y();

   // Check if the curve end has a boundary condition in x.
    if (ps_x1 == ARR_LEFT_BOUNDARY)
    {
      if (ps_x2 == ARR_LEFT_BOUNDARY)
      {
        // Both defined on the left boundary - compare them there.
        CGAL_assertion (ind == ARR_MIN_END);

        return m_traits->compare_y_curve_ends_2_object() (cv, e2->curve(), ind);
      }

      // The curve end is obviously smaller.
      return (SMALLER);
    }

    if (ps_x1 == ARR_RIGHT_BOUNDARY)
    {
      if (ps_x2 == ARR_RIGHT_BOUNDARY)
      {
        // Both defined on the right boundary - compare them there.
        CGAL_assertion (ind == ARR_MAX_END);

        return (m_traits->compare_y_curve_ends_2_object() (cv, e2->curve(), ind));
      }

      // The curve end is obviously larger.
      return (LARGER);
    }

    // ps_x1 == ARR_INTERIOR now!
    CGAL_assertion (ps_x1 == ARR_INTERIOR);

    // Check if the event has a boundary condition in x. Note that if it
    // lies on the left condition, the curve end is larger than it,
    // and if it has a positive boundary condition, the curve end is smaller.
    if (ps_x2 == ARR_LEFT_BOUNDARY)
      return (LARGER);

    if (ps_x2 == ARR_RIGHT_BOUNDARY)
      return (SMALLER);

    // ps_x2 == ARR_INTERIOR now!
    CGAL_assertion (ps_x2 == ARR_INTERIOR);

    // IMPORTANT REMARK HERE: the curve-end is OPEN and ps_x1 == ARR_INTERIOR, thus ps_y1 must be BOTTOM or TOP!
    CGAL_assertion (ps_y1 != ARR_INTERIOR);

    Comparison_result res;

   // Act according to the parameter space of the event.
    if (ps_y2 != ARR_INTERIOR)
    {
      // e2 is an open event here, and thus at most one curve is incident
      Arr_curve_end ind2 = _curve_end(e2);

      // Compare the x-positions of the two entities.
      res = m_traits->compare_x_curve_ends_2_object() (cv, ind,
                                                       e2->curve(), ind2);
      // TODO here e2 can also be a closed event in which case we need to call compare_x_point_curve_end_2

      if (res != EQUAL)
        return (res);

      // In case of equal x-positions, the curve end, compare the boundary conditions on both
      if (ps_y1 == ps_y2) {
        // compare_x_curve_ends for open sides first calls compare_x_on_boundary, and if equal calls compare_x_near_boundary
        return (EQUAL);
      }

      return (ps_y2 == ARR_BOTTOM_BOUNDARY ? LARGER : SMALLER);
    }

    // if we reach here e2 is INTERIOR and thus associated with a valid point.
    CGAL_assertion(ps_y2 == ARR_INTERIOR);
    // while (xcv, ind) still lies on the bottom or top boundary,
    // We compare the x-position of this point with the curve end:
    CGAL_assertion(e2->is_closed());
    res = m_traits->compare_x_point_curve_end_2_object() (e2->point(), cv, ind);

    if (res != EQUAL) {
      // and we need to flip the result
      return (CGAL::opposite(res));
    }

    // In case of equal x-positions, if the curve end lies on the bottom boundary
    // then it lies below e2. Otherwise, it lies on the top boundary above e2.
    return (ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;

  }

  /*! Detemine if the given event represents a left or a right curve end. */
  inline Arr_curve_end _curve_end (const Event* e) const
  {
    // e must be an open event here, and thus exactly one curve is incident
    CGAL_precondition(!e->is_closed());
    return (e->has_left_curves()) ?
      ((e->is_right_end()) ? ARR_MAX_END : ARR_MIN_END) :
      ((e->is_left_end()) ? ARR_MIN_END : ARR_MAX_END);
  }
};

// Forward decleration:
template <class Traits, class Subcurve>
class Sweep_line_event;

/*! \class
 * A functor used to compare curves and curve endpoints by their vertical
 * y-order. Used to maintain the order of the status line (the Y-structure)
 * in the sweep-line algorithm.
 */
template <class Traits_, class Subcurve_>
class Curve_comparer
{
public:

  typedef Traits_                                        Traits_2;
  typedef Subcurve_                                      Subcurve;
  typedef Arr_traits_basic_adaptor_2<Traits_2>           Traits_adaptor_2;

  typedef typename Traits_adaptor_2::Point_2 Point_2;
  typedef typename Traits_adaptor_2::X_monotone_curve_2  X_monotone_curve_2;
  typedef Sweep_line_event<Traits_2, Subcurve>           Event;

private:

  const Traits_adaptor_2 * m_traits;    // A geometric-traits object.
  Event            **m_curr_event;      // Points to the current event point.

public:

  /*! Constructor. */
  template <class Sweep_event>
  Curve_comparer (const Traits_adaptor_2 * t, Sweep_event** e_ptr) :
    m_traits(t),
    m_curr_event(reinterpret_cast<Event**>(e_ptr))
  {}

  /*!
   * Compare the vertical position of two subcurves in the status line.
   * This operator is called only in debug mode.
   */
  Comparison_result operator()(const Subcurve *c1, const Subcurve *c2) const
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
      return (m_traits->compare_y_at_x_right_2_object()
              (c1->last_curve(), c2->last_curve(), (*m_curr_event)->point()));
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
    CGAL_assertion (ps_x1 != ARR_RIGHT_BOUNDARY);

    if (ps_x1 == ARR_LEFT_BOUNDARY)
      return (LARGER);

    // For similar reasons, if c1 begins on the bottom boundary it is below
    // c2, if it is on the top boundary it is above it.
    CGAL_assertion (ps_y1 != ARR_INTERIOR);
    return (ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
  }

  /*!
   * Compare the relative y-order of the given point and the given subcurve.
   */
  Comparison_result operator() (const Point_2& pt, const Subcurve *sc) const
  {
    return (m_traits->compare_y_at_x_2_object()(pt, sc->last_curve()));
  }

};

} //namespace CGAL

#endif
