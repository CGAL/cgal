// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_UTILITY_H
#define CGAL_MCFSKEL_UTILITY_H

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Utility.h
 * @brief This file contains some helper functions like splitting an edge at a
 * given point.
 */

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <cmath>

namespace CGAL {
namespace internal {

template<class TriangleMesh, class TriangleMeshPointPMap, class Traits>
double get_surface_area(TriangleMesh& hg, TriangleMeshPointPMap& hg_point_pmap, const Traits& traits)
{
  typedef typename Traits::Point_3                                       Point;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor    face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  double total_area = 0;
  BOOST_FOREACH(face_descriptor fd, faces(hg))
  {
    halfedge_descriptor hd = halfedge(fd, hg);

    vertex_descriptor v1 = target(hd, hg);
    hd = next(hd, hg);
    vertex_descriptor v2 = target(hd, hg);
    hd = next(hd, hg);
    vertex_descriptor v3 = target(hd, hg);
    Point p1 = get(hg_point_pmap, v1);
    Point p2 = get(hg_point_pmap, v2);
    Point p3 = get(hg_point_pmap, v3);
    total_area += traits.compute_area_3_object()(p1, p2, p3);
  }
  return total_area;
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_UTILITY_H
