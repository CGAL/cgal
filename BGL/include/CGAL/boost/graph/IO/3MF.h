// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé

#ifndef CGAL_BGL_IO_3MF_H
#define CGAL_BGL_IO_3MF_H

#include <CGAL/IO/3MF.h>

#include <string>
#include <vector>

namespace CGAL {

/*!
 * \brief writes the triangle meshes contained in `gs` into the 3mf file `file_name`.
 *
 * \tparam FaceGraphRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose `value type` is
 * a model of the concepts `FaceListGraph` and `HalfedgeListGraph`
 * that has only triangle faces.
 *
 * \param file_name the name of the 3mf file to write.
 * \param gs a `FaceGraphRange` that contains the meshes
 *  to write. An internal property map for `CGAL::vertex_point_t`
 * must be available for each mesh.
 * \param names will contains the name of each mesh in `file_name`.
 *
 * \return `true` if the writing is successful, `false` otherwise.
 */
template<typename FaceGraphRange>
bool write_triangle_meshes_to_3mf(const std::string& file_name,
                                  const FaceGraphRange& gs,
                                  const std::vector<std::string>& names)
{
  typedef typename FaceGraphRange::value_type                                   FaceGraph;
  typedef typename boost::property_map<FaceGraph, boost::vertex_point_t>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  typedef boost::graph_traits<FaceGraph>::vertex_descriptor                     vertex_descriptor;
  typedef boost::graph_traits<FaceGraph>::face_descriptor                       face_descriptor;

  typedef std::vector<std::size_t>                                              Polygon;
  typedef std::vector<Polygon>                                                  PolygonRange;
  typedef std::vector<Point>                                                    PointRange;

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;

  for(const auto& g : gs)
  {
    PointRange points;
    points.reserve(num_vertices(g));
    PolygonRange triangles;
    triangles.reserve(num_faces(g));

    VPM vpm = get(boost::vertex_point, g);
    std::unordered_map<typename boost::graph_traits<FaceGraph>::vertex_descriptor, std::size_t> vertex_id_map;

    std::size_t i = 0;
    for(const vertex_descriptor v : vertices(g))
    {
      points.push_back(get(vpm, v));
      vertex_id_map[v] = i++;
    }

    all_points.push_back(points);
    for(const face_descriptor f : faces(g))
    {
      Polygon triangle;
      for(vertex_descriptor vert : CGAL::vertices_around_face(halfedge(f, g), g))
        triangle.push_back(vertex_id_map[vert]);

      triangles.push_back(triangle);
    }

    all_polygons.push_back(triangles);
  }

  return write_triangle_soups_to_3mf(file_name, all_points, all_polygons, names);
}

} // namespace CGAL

#endif // CGAL_BGL_IO_3MF_H
