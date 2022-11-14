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
#include <CGAL/Arr_tags.h>
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

  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;

  /*!
   */
  typedef typename Arr_sides_category<Left_side_category,
                                      Bottom_side_category,
                                      Top_side_category,
                                      Right_side_category>::result
    Sides_category;

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

    if ((ps_x1 == ARR_INTERIOR) && (ps_y1 == ARR_INTERIOR)) {
      // The first curve has a valid left endpoint. Compare the y-position
      // of this endpoint with the second subcurve.
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
    // c2, if it begins on the top boundary it is above it.
    CGAL_assertion (ps_y1 != ARR_INTERIOR);
    return (ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
  }

  /*! Compare the relative y-order of the given point and the given subcurve.
   */
  Comparison_result operator()(const Point_2& pt, const Subcurve* sc) const
  { return compare_point_curve(pt, sc->last_curve(), Sides_category()); }

private:
  Comparison_result compare_point_curve(const Point_2& pt,
                                        const X_monotone_curve_2& cv,
                                        Arr_all_sides_oblivious_tag) const
  { return m_traits->compare_y_at_x_2_object()(pt, cv); }

  Comparison_result compare_point_curve(const Point_2& pt,
                                        const X_monotone_curve_2& cv,
                                        Arr_all_sides_not_finite_tag) const
  { return m_traits->compare_y_at_x_2_object()(pt, cv); }

  Comparison_result compare_point_curve(const Point_2& pt,
                                        const X_monotone_curve_2& cv,
                                        Arr_not_all_sides_not_finite_tag) const
  {
    Arr_parameter_space ps_y1 = (*m_curr_event)->parameter_space_in_y();

#if 0
    Arr_parameter_space ps_x1 = (*m_curr_event)->parameter_space_in_x();
    typename Traits_adaptor_2::Parameter_space_in_x_2 ps_x_op =
      m_traits->parameter_space_in_x_2_object();
    Arr_parameter_space ps_x_min = ps_x_op(cv, ARR_MIN_END);
    Arr_parameter_space ps_x_max = ps_x_op(cv, ARR_MAX_END);
    typename Traits_adaptor_2::Parameter_space_in_y_2 ps_y_op =
      m_traits->parameter_space_in_y_2_object();
    Arr_parameter_space ps_y_min = ps_y_op(cv, ARR_MIN_END);
    Arr_parameter_space ps_y_max = ps_y_op(cv, ARR_MAX_END);

    CGAL::set_pretty_mode(std::cout);
    std::cout << "\n FUNCTOR pt-cv" << std::endl;
    std::cout << "pt: " << pt << std::endl;
    std::cout << "ps_x1: " << ps_x1 << std::endl;
    std::cout << "ps_y1: " << ps_y1 << std::endl;
    std::cout << "cv: " << cv << std::endl;
    std::cout << "ps_x_min: " << ps_x_min << std::endl;
    std::cout << "ps_y_min: " << ps_y_min << std::endl;
    std::cout << "ps_x_max: " << ps_x_max << std::endl;
    std::cout << "ps_y_max: " << ps_y_max << std::endl;
#endif

    if (ps_y1 == ARR_TOP_BOUNDARY) return LARGER;
    if (ps_y1 == ARR_BOTTOM_BOUNDARY) return SMALLER;
    return m_traits->compare_y_at_x_2_object()(pt, cv);
  }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
