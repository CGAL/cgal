// Copyright (c) 2013 INRIA Sophia-Anitpolis (France).
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
//
//
// Author(s)     : Pierre Alliez


#ifndef CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
#define CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H

#include <CGAL/boost/graph/helpers.h>

namespace CGAL{

namespace Polygon_mesh_processing{

/**
* \ingroup PkgPolygonMeshProcessing
* computes the outward unit vector normal to facet `f`.
* @tparam Kernel a %CGAL `Kernel` with `FT` a model of `FieldWithSqrt`
* @tparam PolygonMesh a model of `FaceGraph`
*
* @param f the facet on which the normal is computed
* @param pmesh the polygon mesh to which `f` belongs
*/
template <class Kernel, class PolygonMesh>
typename Kernel::Vector_3
compute_facet_normal(
  typename boost::graph_traits<PolygonMesh>::face_descriptor f,
  const PolygonMesh& pmesh)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type
    ppmap = get(boost::vertex_point, pmesh);

  Vector normal = CGAL::NULL_VECTOR;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  halfedge_descriptor he = halfedge(f, pmesh);
  halfedge_descriptor end = he;
  do
  {
    const Point& prv = ppmap[target(prev(he, pmesh), pmesh)];
    const Point& curr = ppmap[target(he, pmesh)];
    const Point& nxt = ppmap[target(next(he, pmesh), pmesh)];
    Vector n = CGAL::cross_product(nxt - curr, prv - curr);
    normal = normal + n;

    he = next(he, pmesh);
  } while (he != end);

  return normal / std::sqrt(normal * normal);
}

/**
* \ingroup PkgPolygonMeshProcessing
* computes the unit normal at vertex `v` as the average of the normals of incident facets.
* @tparam Kernel a %CGAL `Kernel` with `FT` a model of `FieldWithSqrt`
* @tparam PolygonMesh a model of `FaceGraph`
*
* @param v the vertex around which the normal is computed
* @param pmesh the polygon mesh to which `v` belongs
*/
template<typename Kernel, typename PolygonMesh>
typename Kernel::Vector_3
compute_vertex_normal(
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor v,
    const PolygonMesh& pmesh)
{
  typedef typename Kernel::Vector_3 Vector;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  Vector normal = CGAL::NULL_VECTOR;
  halfedge_descriptor he = halfedge(v, pmesh);
  halfedge_descriptor end = he;
  do
  {
    if (!is_border(he, pmesh))
    {
      Vector n = compute_facet_normal<Kernel>(face(he, pmesh), pmesh);
      normal = normal + (n / std::sqrt(n*n));
    }
    he = opposite(next(he, pmesh), pmesh);
  } while (he != end);

  return normal / std::sqrt(normal * normal);
}

} } // end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_COMPUTE_NORMAL_H
