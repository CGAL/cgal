// Copyright (c) 2014 GeometryFactory
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
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_MISC_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_MISC_H

#include <CGAL/license/Surface_mesh_shortest_path.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

namespace internal {

template <class Triangle_3, class Triangle_mesh, class VertexPointMap>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor edge, const Triangle_mesh& g, VertexPointMap vertexPointMap)
{
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor e0 = edge;
  halfedge_descriptor e1 = next(edge, g);

  return Triangle_3(get(vertexPointMap, source(e0, g)), get(vertexPointMap, target(e0, g)), get(vertexPointMap, target(e1, g)));
}

template <class Triangle_3, class Triangle_mesh>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor edge, const Triangle_mesh& g)
{
  return triangle_from_halfedge<Triangle_3, Triangle_mesh, typename boost::property_map<Triangle_mesh, boost::vertex_point_t>::type>(edge, g, get(boost::vertex_point, g));
}


template <class Triangle_mesh>
size_t edge_index(typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor he, Triangle_mesh& p)
{
  typedef typename boost::graph_traits<Triangle_mesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  face_descriptor f = face(he, p);

  halfedge_descriptor start = halfedge(f, p);
  halfedge_descriptor current = start;

  size_t count = 0;

  while (current != he)
  {
    current = next(current, p);
    ++count;
  }

  return count;
}

template <class FT>
FT internal_sqrt(const FT& x, CGAL::Tag_true)
{
  return CGAL::sqrt(x);
}

template <class FT>
FT internal_sqrt(const FT& x, CGAL::Tag_false)
{
  return FT(std::sqrt(CGAL::to_double(x)));
}

template <class FT>
FT select_sqrt(const FT& x)
{
  typedef ::CGAL::Algebraic_structure_traits<FT> AST;
  static const bool has_sqrt = ! ::boost::is_same< ::CGAL::Null_functor, typename AST::Sqrt >::value;
  return internal_sqrt(x, ::CGAL::Boolean_tag<has_sqrt>());
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_MISC_H
