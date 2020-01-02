// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/intersections.h>

#include <boost/optional.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {

template<class Placement, class Tree>
class Bounded_distance_placement
{
public:
  typedef typename Placement::TM                          TM;

public:
  Bounded_distance_placement(const double dist,
                             const Tree& tree,
                             const Placement& placement = Placement())
    : m_get_placement(placement), m_tree(tree), m_threshold_dist(dist)
  { }

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()(const Profile& profile) const
  {
    typedef typename Profile::Point                       Point;
    typedef typename CGAL::Kernel_traits< Point >::Kernel Kernel;

    boost::optional<typename Profile::Point> op = m_get_placement(profile);
    if(op)
    {
      CGAL_assertion(!m_tree.empty());

      const Point& p = *op;

      m_tree.accelerate_distance_queries();
      const Point& cp = m_tree.best_hint(p).first; // requires accelerate distance query to be called.

      // We could do better by having access to the internal kd-tree
      // and call search_any_point with a fuzzy_sphere.
      const double sqtd = CGAL::square(m_threshold_dist);

      // if no input vertex is closer than the threshold, then
      // any face closer than the threshold is intersected by
      // the sphere (avoid the inclusion of the mesh into the threshold sphere)
      if(CGAL::compare_squared_distance(p, cp, sqtd) != LARGER ||
         m_tree.do_intersect(CGAL::Sphere_3<Kernel>(p, sqtd)))
        return op;

      return boost::optional<Point>();
    }

    return op;
  }

private:
  const Placement m_get_placement;
  const Tree& m_tree;
  const double m_threshold_dist;
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H
