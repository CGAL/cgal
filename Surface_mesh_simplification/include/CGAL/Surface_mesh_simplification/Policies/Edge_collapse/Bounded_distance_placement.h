// Copyright (c) 2017  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include <boost/optional.hpp>
#include<CGAL/intersections.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

template<class Placement>
class Bounded_distance_placement
{
public:

  typedef typename Placement::TM TM ;

public:

  Bounded_distance_placement(const double sq_dist,
                             const Placement& placement = Placement() )
    : mPlacement(placement), threshold_sq_dist(sq_dist)
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
      const Point* p = boost::get<Point>(&op);
      const typename Profile::Triangle_vector& triangles = aProfile.triangles();
      typename Profile::Triangle_vector::const_iterator it = triangles.begin();
      typename Profile::VertexPointMap ppmap = aProfile.vertex_point_map();

      bool does_intersect = false;
      CGAL::Sphere_3<Kernel> s(*p, threshold_sq_dist);
      while(it!= triangles.end()){
        typename Kernel::Triangle_3 t(get(ppmap, it->v0), get(ppmap, it->v1), get(ppmap, it->v2));
        if(CGAL::do_intersect(t, s)){
          does_intersect = true;
          break;
        }
        ++it;
      }
      if(!does_intersect)
      {
        return boost::none;
      }
    }
    return op;
  }


private:

  Placement  mPlacement ;
  double threshold_sq_dist;

};
}
}

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_DISTANCE_PLACEMENT_H


