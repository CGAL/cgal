// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Maxime Gimeno
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include <boost/optional.hpp>
#include<CGAL/intersections.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

template<class Placement, class Tree>
class Bounded_distance_placement
{
public:

  typedef typename Placement::TM TM ;

public:

  Bounded_distance_placement(const double dist,
                              const Tree& tree,
                             const Placement& placement = Placement() )
    : mPlacement(placement), tree(tree), threshold_dist(dist)
  {
  }

  template <typename Profile>
  boost::optional<typename Profile::Point>
  operator()( Profile const& aProfile) const
  {
    boost::optional<typename Profile::Point> op = mPlacement(aProfile);
    typedef typename Profile::Point Point;
    typedef typename CGAL::Kernel_traits< Point >::Kernel Kernel;
    if(op){

      CGAL_assertion(!tree.empty());

      const Point* p = boost::get<Point>(&op);

      tree.accelerate_distance_queries();
      Point cp = tree.best_hint(*p).first; // requires accelerate distance query to be called.
                                    // We could do better by having access to the internal kd-tree
                                    // and call search_any_point with a fuzzy_sphere.

      const double sqtd = threshold_dist*threshold_dist;

      // if no input vertex is closer than the threshold, then
      // any face closer than the threshold is intersected by
      // the sphere (avoid the inclusion of the mesh into the threshold sphere)
      if( CGAL::compare_squared_distance(*p,cp, sqtd)!=LARGER ||
          tree.do_intersect(CGAL::Sphere_3<Kernel>(*p, sqtd)))
      {
        return op;
      }
      return boost::none;
    }
    return op;
  }


private:

  Placement  mPlacement ;
  const Tree& tree;
  double threshold_dist;

};
}
}

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H


