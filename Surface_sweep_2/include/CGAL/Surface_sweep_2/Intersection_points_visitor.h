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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_INTERSECTION_POINTS_VISITOR_H
#define CGAL_SURFACE_SWEEP_2_INTERSECTION_POINTS_VISITOR_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of a surface-sweep visitor that reports all intersection points
 * among a set of input curves.
 */

#include <vector>
// #include <iterator>

#include <CGAL/Surface_sweep_2/Default_visitor.h>
#include <CGAL/Surface_sweep_2/Surface_sweep_2_utils.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Intersection_points_visitor
 *
 * A simple surface-sweep visitor that reports all intersection points among a
 * set of input curves. Used by compute_intersection_points().
 */
template <typename GeometryTraits_2, typename OutputIterator,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Intersection_points_visitor :
  public Default_visitor<Intersection_points_visitor<GeometryTraits_2,
                                                      OutputIterator,
                                                     Allocator_>,
                         GeometryTraits_2, Allocator_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef OutputIterator                                Output_iterator;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Intersection_points_visitor<Gt2, Output_iterator, Allocator>
                                                        Self;
  typedef Default_visitor<Self, Gt2, Allocator>         Base;

public:
  typedef typename Base::Event                          Event;
  typedef typename Base::Subcurve                       Subcurve;

  typedef typename Subcurve::Status_line_iterator       Status_line_iterator;

  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Base::Surface_sweep_2                Surface_sweep_2;

protected:
  Output_iterator m_out;                 // The output points.
  bool m_includeEndPoints;              // Should we include endpoints.

public:
  Intersection_points_visitor(Output_iterator out, bool endpoints) :
    m_out(out),
    m_includeEndPoints(endpoints)
  {}

  template <typename CurveIterator>
  void sweep(CurveIterator begin, CurveIterator end)
  {
    std::vector<X_monotone_curve_2> curves_vec;
    std::vector<Point_2> points_vec;

    curves_vec.reserve(std::distance(begin,end));
    make_x_monotone(begin, end,
                    std::back_inserter(curves_vec),
                    std::back_inserter(points_vec),
                    this->traits());

    //Perform the sweep
    Surface_sweep_2* sl = this->surface_sweep();
    sl->sweep(curves_vec.begin(), curves_vec.end(),
              points_vec.begin(), points_vec.end());
  }

  bool after_handle_event(Event* event,
                          Status_line_iterator /* iter */,
                          bool /* flag */)
  {
    if ((m_includeEndPoints ||
         event->is_intersection() ||
         event->is_weak_intersection()) && event->is_closed())
    {
      *m_out = event->point();
      ++m_out;
    }
    return true;
  }

  Output_iterator output_iterator() { return m_out; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
