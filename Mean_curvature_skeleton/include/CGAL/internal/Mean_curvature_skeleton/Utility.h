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

/**
* Split the edge
* @param hg the mesh containing the given edge
* @param ei the edge to be split
* @param pn the position of the new vertex created by the split
*/
template<class TriangleMesh, class TriangleMeshPointPMap>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
mesh_split(TriangleMesh& hg,
           TriangleMeshPointPMap& hg_point_pmap,
           typename boost::graph_traits<TriangleMesh>::halfedge_descriptor ei,
           typename TriangleMesh::Traits::Point_3 pn)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor            halfedge_descriptor;


  // halfedge_descriptor en = Euler::split_edge(ei, hg); // there is an issue in this function for now use the polyhedron version in the meantime
  halfedge_descriptor en = hg.split_edge(ei);
  en->vertex()->vertices.clear();
  boost::put(hg_point_pmap, target(en,hg), pn);
  Euler::split_face(en, next(ei,hg), hg);

  en->id() = -1;
  opposite(en,hg)->id() = -1;
  ei->id() = -1;
  opposite(ei,hg)->id() = -1;
  next(en,hg)->id() = -1;
  opposite(next(en,hg),hg)->id() = -1;
  next(next(en,hg),hg)->id() = -1;
  next(en,hg)->id() = -1; // AF: the same as 3 lines above? error or duplicate?
  halfedge_descriptor ej = opposite(en, hg);
  if (! is_border(ej,hg))
  {
    Euler::split_face(opposite(ei,hg), next(ej,hg), hg);
    next(ej,hg)->id() = -1;
    halfedge_descriptor ei_op_next = next(opposite(ei,hg),hg);
    ei_op_next->id() = -1;
    opposite(ei_op_next,hg)->id() = -1;
    next(ei_op_next,hg)->id() = -1;
  }

  return en;
}

template<class Vertex, class Kernel>
double get_triangle_area(typename Kernel::Point_3 p1,
                         typename Kernel::Point_3 p2,
                         typename Kernel::Point_3 p3)
{
  typedef typename Kernel::Vector_3 Vector;
  Vector v12(p1, p2);
  Vector v13(p1, p3);
  return std::sqrt(cross_product(v12, v13).squared_length()) * 0.5;
}

template<class TriangleMesh, class TriangleMeshPointPMap>
double get_surface_area(TriangleMesh& hg, TriangleMeshPointPMap& hg_point_pmap)
{
  typedef typename TriangleMesh::Traits                                  Kernel;
  typedef typename Kernel::Point_3                                       Point;
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
    Point p1 = boost::get(hg_point_pmap, v1);
    Point p2 = boost::get(hg_point_pmap, v2);
    Point p3 = boost::get(hg_point_pmap, v3);
    total_area += internal::get_triangle_area<vertex_descriptor, Kernel>(p1, p2, p3);
  }
  return total_area;
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_UTILITY_H
