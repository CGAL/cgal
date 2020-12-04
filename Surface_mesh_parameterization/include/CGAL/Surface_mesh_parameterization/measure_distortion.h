// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_DISTORTION_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_DISTORTION_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/circulator.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <vector>

namespace CGAL {
namespace Surface_mesh_parameterization {

#ifndef DOXYGEN_RUNNING

// Measure L2 stretch
template <typename VertexRange, typename FaceRange, typename TriangleMesh, typename VertexUVmap>
double compute_L2_stretch(const VertexRange& vertex_range,
                          const FaceRange& face_range,
                          const TriangleMesh& tmesh,
                          const VertexUVmap uvmap)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type             Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::Kernel                           Kernel;
  typedef typename Kernel::Point_2                                                Point_2;

  typedef CGAL::dynamic_face_property_t<double>                                   Face_double_tag;
  typedef typename boost::property_map<TriangleMesh, Face_double_tag>::const_type Face_double_map;

  Face_double_map area_2D = get(Face_double_tag(), tmesh);
  Face_double_map area_3D = get(Face_double_tag(), tmesh);

  // iterate fpr all inner vertices and for each vertex
  std::vector<double> area_dist;

  double A_3D = 0.;
  double A_2D = 0.;

  for(face_descriptor f : face_range)
  {
    std::vector<Point_2> uv_points;
    for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh), tmesh))
      uv_points.push_back(get(uvmap, v));

    const double a_2D = abs(CGAL::area(get(uvmap, target(halfedge(f, tmesh), tmesh)),
                                       get(uvmap, source(halfedge(f, tmesh), tmesh)),
                                       get(uvmap, target(next(halfedge(f, tmesh), tmesh), tmesh))));
    const double a_3D = Polygon_mesh_processing::face_area(f, tmesh);

    put(area_2D, f, a_2D);
    put(area_3D, f, a_3D);

    A_2D += a_2D;
    A_3D += a_3D;
  }

  for(vertex_descriptor v : vertex_range)
  {
    // inner vertices only
    if(CGAL::is_border(v, tmesh))
      continue;

    double a_2D = 0.;
    double a_3D = 0.;

    // find the area of all the adjacent faces to this vertex
    CGAL::Face_around_target_circulator<TriangleMesh> f_j(halfedge(v, tmesh), tmesh), end = f_j;
    CGAL_For_all(f_j, end)
    {
      if(*f_j == boost::graph_traits<TriangleMesh>::null_face())
        continue;

      a_2D += get(area_2D, *f_j);
      a_3D += get(area_3D, *f_j);
    }

    a_2D /= A_2D;
    a_3D /= A_3D;

    area_dist.push_back(square((a_3D/a_2D) - 1.));
  }

  return sqrt(std::accumulate(area_dist.begin(), area_dist.end(), 0.));
}

template <typename TriangleMesh, typename VertexUVmap>
double compute_L2_stretch(const TriangleMesh& tmesh,
                          const VertexUVmap uvmap)
{
  return compute_L2_stretch(vertices(tmesh), faces(tmesh), tmesh, uvmap);
}

#endif // DOXYGEN_RUNNING

} // namespace Surface_mesh_parameterization
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_DISTORTION_H
