// Copyright (c) 2021 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIM_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIM_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

/**
  * \ingroup PMP_corefinement_grp
  *
  * trims `tm` at a polyline that must be close to the surface mesh.
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph`, and `FaceListGraph`.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  * @tparam PointRange a range of 3D points
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tm input triangulated surface mesh that will also contain one output mesh
  * @param other the second output mesh
  * @param pr the point range for a closed polyline  that is close to `tm`. The first and last point must be equal.
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
**/

template <class TriangleMesh,
          class PointRange,
          class NamedParameters = parameters::Default_named_parameters>
bool trim(TriangleMesh& tm,
          TriangleMesh& other,
          const PointRange& pr,
          const NamedParameters& np = parameters::default_values())
{

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type     Geom_traits;
  typedef typename Geom_traits::Point_3                                   Point_3;
  typedef typename Geom_traits::Vector_3                                  Vector_3;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh>          AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<Geom_traits, AABB_face_graph_primitive>       AABB_face_graph_traits;
  typedef CGAL::AABB_tree<AABB_face_graph_traits>                         Tree;
  typedef PMP::Face_location<TriangleMesh, double>                        Face_location;

  Tree tree;
  build_AABB_tree(tm, tree);

  int n = std::distance(pr.begin(), pr.end());
  const double extrusion_distance = 0.1;
  typename PointRange::const_iterator it = pr.begin();
  std::vector<Point_3> points;
  points.reserve(2* (n-1));

  for(int i = 0; i < n-1; ++i){
    const Point_3 p = *it;
    ++it;
    Face_location location = locate_with_AABB_tree(p, tree, tm);
    face_descriptor fd = location.first;
    Vector_3 v = compute_face_normal(fd, tm);
    v /= sqrt(v.squared_length());
    v *= extrusion_distance;
    Point_3 q = construct_point(location, tm);
    points.push_back(q+v);
    points.push_back(q-v);
  }
  std::vector<std::array<int,3>> faces;
  for(int i = 0; i < (2*(n-1)-2); ++i){
    std::array<int,3> a = {i, i+1, i+2};
    faces.push_back(a);
  }

  // We could have used a modulo for closing
  {
    std::array<int,3> a = {2*(n-1)-2, 2*(n-1)-1, 0};
    faces.push_back(a);
  }
  {
    std::array<int,3> a = {2*(n-1)-1, 0, 1};
    faces.push_back(a);
  }

  // The extruded polyline may have problems
  bool orientable = orient_polygon_soup(points,faces);
  if(! orientable){
    std::cerr << "Provide better curve" << std::endl;
    return false;
  }

  Mesh splitter;
  polygon_soup_to_polygon_mesh(points, faces, splitter);

  // CGAL::IO::write_polygon_mesh("splitter.off", splitter);

  split(tm, splitter);

  std::vector<Mesh> cctm;

  split_connected_components(tm, cctm);

  if(cctm.size() != 2){
    std::cerr << "We expected 2 components and not " << cctm.size() << std::endl;
    return false;
  }

  tm = cctm[0];
  other = cctm[1];
  return true;
}


} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIM_H
