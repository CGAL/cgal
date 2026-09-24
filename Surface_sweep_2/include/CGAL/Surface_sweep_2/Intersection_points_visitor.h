// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Ron Wein <wein@post.tau.ac.il>
//             Efi Fogel <efif@post.tau.ac.il>

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
  public Default_visitor<Intersection_points_visitor<GeometryTraits_2, OutputIterator, Allocator_>,
                         GeometryTraits_2, Allocator_> {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Output_iterator = OutputIterator;
  using Allocator = Allocator_;

private:
  using Gt2 = Geometry_traits_2;
  using Self = Intersection_points_visitor<Gt2, Output_iterator, Allocator>;
  using Base = Default_visitor<Self, Gt2, Allocator>;

public:
  using Event = typename Base::Event;
  using Subcurve = typename Base::Subcurve;

  using Status_line_iterator = typename Subcurve::Status_line_iterator;

  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;
  using Point_2 = typename Gt2::Point_2;

  using Surface_sweep_2 = typename Base::Surface_sweep_2;

protected:
  Output_iterator m_out;                // the output points.
  bool m_includeEndPoints;              // should we include endpoints.

public:
  /*! constructs.
   */
  Intersection_points_visitor(Output_iterator out, bool endpoints) :
    m_out(out),
    m_includeEndPoints(endpoints)
  {}

  /*!
   */
  bool after_handle_event(Event* event, Status_line_iterator /* iter */, bool /* flag */) {
    if ((m_includeEndPoints || event->is_intersection() || event->is_weak_intersection()) && event->is_closed()) {
      *m_out = event->point();
      ++m_out;
    }
    return true;
  }

  /*!
   */
  Output_iterator output_iterator() { return m_out; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
