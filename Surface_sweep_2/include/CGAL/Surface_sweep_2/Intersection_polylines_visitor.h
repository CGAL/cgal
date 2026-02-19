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

#ifndef CGAL_SURFACE_SWEEP_2_INTERSECTION_POLYLINES_VISITOR_H
#define CGAL_SURFACE_SWEEP_2_INTERSECTION_POLYLINES_VISITOR_H

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

/*! \class Intersection_polylines_visitor
 *
 * A simple surface-sweep visitor that reports all intersection points among a
 * set of input curves. Used by compute_intersection_points().
 */
template <typename GeometryTraits_2, typename Points_, typename Polylines_, typename Allocator_ = CGAL_ALLOCATOR(int)>
class Intersection_polylines_visitor :
    public Default_visitor<Intersection_polylines_visitor<GeometryTraits_2, Points_, Polylines_, Allocator_>,
                           GeometryTraits_2, Allocator_> {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Points = Points_;
  using Polylines = Polylines_;
  using Allocator = Allocator_;

private:
  using Gt2 = Geometry_traits_2;
  using Self = Intersection_polylines_visitor<Gt2, Points, Polylines, Allocator>;
  using Base = Default_visitor<Self, Gt2, Allocator>;

public:
  using Event = typename Base::Event;
  using Subcurve = typename Base::Subcurve;
  using Status_line_iterator = typename Subcurve::Status_line_iterator;
  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;
  using Point_2 = typename Gt2::Point_2;
  using Surface_sweep_2 = typename Base::Surface_sweep_2;

protected:
  Points& m_points;
  Polylines& m_polylines;

public:
  Intersection_polylines_visitor(Points& points, Polylines& polylines) :
    m_points(points),
    m_polylines(polylines)
  {}

  /*!
   */
  template <typename CurveIterator>
  void sweep(CurveIterator begin, CurveIterator end) {
    m_points.reserve(2 * std::distance(begin,end));
    m_polylines.resize(std::distance(begin,end));
    Surface_sweep_2* sl = this->surface_sweep();
    sl->sweep(begin, end);
  }

  /*!
   */
  bool after_handle_event(Event* event, Status_line_iterator /* iter */, bool /* flag */) {
    if (! event->is_closed()) return true;
    if (! event->has_left_curves() && ! event->has_right_curves()) return true;

    auto id = m_points.size();
    m_points.push_back(event->point());

    if (! event->has_right_curves() && event->has_left_curves()) {
      for (auto it = event->left_curves_begin(); it != event->left_curves_end(); ++it) {
        const auto* sc = *it;
        m_polylines[sc->last_curve().data()].push_back(id);
      }
      return true;
    }

    if (! event->has_left_curves() && event->has_right_curves()) {
      for (auto it = event->right_curves_begin(); it != event->right_curves_end(); ++it) {
        const auto* sc = *it;
        m_polylines[sc->last_curve().data()].push_back(id);
      }
      return true;
    }

    std::unordered_set<std::size_t> cids;
    for (auto it = event->left_curves_begin(); it != event->left_curves_end(); ++it) {
      const auto* sc = *it;
      cids.insert(sc->last_curve().data());
    }
    for (auto it = event->right_curves_begin(); it != event->right_curves_end(); ++it) {
      const auto* sc = *it;
      cids.insert(sc->last_curve().data());
    }
    for (auto cid : cids) m_polylines[cid].push_back(id);
    cids.clear();

    return true;
  }

  const Points& points() const { return m_points; }
  Points& points() { return m_points; }

  const Polylines& polylines() const { return m_polylines; }
  Polylines& polylines() { return m_polylines; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
