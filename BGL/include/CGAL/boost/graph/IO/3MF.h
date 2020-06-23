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

#include <CGAL/boost/graph/iterator.h>

#include <boost/property_map/property_map.hpp>

#include <string>
#include <unordered_map>
#include <vector>

#if defined(CGAL_LINKED_WITH_3MF) || defined(DOXYGEN_RUNNING)

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgBGLIOFuncs3MF
 *
 * \brief Writes the triangle meshes contained in `gs` into the 3mf file `filename`.
 *
 * \tparam GraphRange a model of the concepts `RandomAccessContainer`
 *                    and `BackInsertionSequence` whose `value type` is
 *                    a model of the concepts `FaceGraph` and `HalfedgeListGraph`
 *                    that has only triangle faces.
 *
 * \param filename the name of the 3mf file to write.
 * \param gs a container of triangle meshes to write. An internal property map for `CGAL::vertex_point_t`
 *           must be available for each mesh.
 * \param names a range of `std::string` associating a name to each mesh to be written out, which
 *              will appear in the output.
 *
 * \return `true` if the writing is successful, `false` otherwise.
 *
 * \sa `read_3MF()`
 */
template<typename GraphRange>
bool write_3MF(const std::string& filename,
               const GraphRange& gs,
               const std::vector<std::string>& names)
{
  typedef typename GraphRange::value_type                                       FaceGraph;
  typedef typename boost::property_map<FaceGraph, boost::vertex_point_t>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor            vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor              face_descriptor;

  typedef std::vector<std::size_t>                                              Polygon;
  typedef std::vector<Polygon>                                                  PolygonRange;
  typedef std::vector<Point>                                                    PointRange;

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;

  for(const FaceGraph& g : gs)
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

  return write_triangle_soups_to_3mf(filename, all_points, all_polygons, names);
}

} // namespace CGAL

#endif // defined(CGAL_LINKED_WITH_3MF) || defined(DOXYGEN_RUNNING)

#endif // CGAL_BGL_IO_3MF_H
