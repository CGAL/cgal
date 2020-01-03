// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_MISC_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_MISC_H

#include <CGAL/license/Surface_mesh_shortest_path.h>


#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {
namespace Surface_mesh_shortest_paths_3 {
namespace internal {

template <class Triangle_3, class Triangle_mesh, class VertexPointMap>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor edge,
                                  const Triangle_mesh& g,
                                  const VertexPointMap vertexPointMap)
{
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;

  const halfedge_descriptor e0 = edge;
  const halfedge_descriptor e1 = next(edge, g);

  return Triangle_3(get(vertexPointMap, source(e0, g)),
                    get(vertexPointMap, target(e0, g)),
                    get(vertexPointMap, target(e1, g)));
}

template <class Triangle_3, class Triangle_mesh>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor edge,
                                  const Triangle_mesh& g)
{
  return triangle_from_halfedge<Triangle_3>(edge, g, get(boost::vertex_point, g));
}

template <class Triangle_mesh>
std::size_t edge_index(typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor he,
                       const Triangle_mesh& p)
{
  typedef typename boost::graph_traits<Triangle_mesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  const face_descriptor f = face(he, p);

  const halfedge_descriptor start = halfedge(f, p);
  halfedge_descriptor current = start;

  std::size_t count = 0;
  while (current != he)
  {
    current = next(current, p);
    ++count;
  }

  return count;
}

} // namespace internal
} // namespace Surface_mesh_shortest_paths_3
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_MISC_H
