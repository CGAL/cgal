// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_DO_INTERIOR_INTERSECT_VISITORS_H
#define CGAL_SURFACE_SWEEP_2_DO_INTERIOR_INTERSECT_VISITORS_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the basic sweep-line visitors, for the usage of the global
 * sweep-line functions.
 */

#include <vector>

#include <CGAL/Surface_sweep_2/Default_visitor.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Do_interior_intersect_visitor
 *
 * A simple sweep-line visitor that determines whether the curves in a given set
 * intersect in their interiors.
 */
template <typename GeometryTraits_2,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Do_interior_intersect_visitor :
  public Default_visitor<Do_interior_intersect_visitor<GeometryTraits_2,
                                                       Allocator_>,
                         GeometryTraits_2, Allocator_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Do_interior_intersect_visitor<Gt2, Allocator> Self;
  typedef Default_visitor<Self, Gt2, Allocator>         Base;

public:
  typedef typename Base::Event                          Event;
  typedef typename Base::Subcurve                       Subcurve;

  typedef typename Subcurve::Status_line_iterator       Status_line_iterator;

  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Base::Surface_sweep_2                Surface_sweep_2;

protected:
  // Data members:
  bool m_found_x;               // Have we found an intersection so far.

public:
  Do_interior_intersect_visitor() : m_found_x(false) {}

  template <typename CurveIterator>
  void sweep(CurveIterator begin, CurveIterator end)
  {
    std::vector<X_monotone_curve_2> curves_vec;
    std::vector<Point_2> points_vec;

    curves_vec.reserve(std::distance(begin,end));
    make_x_monotone(begin, end,
                    std::back_inserter(curves_vec),
                    std::back_inserter(points_vec),
                    this-> traits());

    // Perform the sweep.
    Surface_sweep_2* sl = this->surface_sweep();
    sl->sweep(curves_vec.begin(), curves_vec.end(),
              points_vec.begin(), points_vec.end());
  }

  void update_event(Event* /* e */,
                    Subcurve* /* sc1 */,
                    Subcurve* /* sc2 */,
                    bool /* is_new */)
  { m_found_x = true; }

  void update_event(Event* /* e */,
                    Subcurve* /* sc1 */)
  { m_found_x = true; }

  void update_event(Event* /* e */,
                    const Point_2& /* end_point */,
                    const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */,
                    bool /* is_new */)
  {}

  void update_event(Event* /* e */,
                    const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */,
                    bool /* is_new */)
  {}

  void update_event(Event* /* e */,
                    const Point_2& /* pt */,
                    bool /* is_new */)
  {}

  template <typename XCurveIterator>
  void sweep_xcurves(XCurveIterator begin, XCurveIterator end)
  {
    // Perform the sweep.
    Surface_sweep_2* sl = this->surface_sweep();
    sl->sweep(begin, end);
  }

  void found_overlap(Subcurve* /* sc1 */,
                     Subcurve* /* sc2 */,
                     Subcurve* /* ov_sc */)
  { m_found_x = true; }

  bool after_handle_event(Event* /* event */,
                          Status_line_iterator /* iter */,
                          bool /* flag */)
  {
    if (m_found_x) {
      Surface_sweep_2* sl = this->surface_sweep();
      sl->stop_sweep();
    }
    return true;
  }

  bool found_intersection() { return m_found_x; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
