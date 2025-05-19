// Copyright (c) 2015-2020  GeometryFactory (France).  All rights reserved.
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
#include <CGAL/IO/helpers.h>

#include <CGAL/boost/graph/iterator.h>

#include <boost/property_map/property_map.hpp>

#include <string>
#include <unordered_map>
#include <vector>
#include <type_traits>

#if defined(CGAL_LINKED_WITH_3MF) || defined(DOXYGEN_RUNNING)

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

/*!
 * \ingroup PkgBGLIoFuncs3MF
 *
 * \brief writes the triangle meshes contained in `gs` into the file `filename`, using the \ref IOStream3MF.
 *
 * \tparam GraphRange a model of the concepts `RandomAccessContainer`
 *                    and `BackInsertionSequence` whose `value_type` is
 *                    a model of the concepts `FaceGraph` and `HalfedgeListGraph`
 *                    that has only triangle faces.
 *
 * \param filename the name of the 3mf file to write
 * \param gs a container of triangle meshes to write. An internal property map for `CGAL::vertex_point_t`
 *           must be available for each mesh.
 * \param names a range of `std::string` associating a name to each mesh to be written out, which
 *              will appear in the output
 *
 * \return `true` if the writing is successful, `false` otherwise.
 *
 * \sa `read_3MF()`
 */
template<typename GraphRange>
bool write_3MF(const std::string& filename,
               const GraphRange& gs,
               const std::vector<std::string>& names
#ifndef DOXYGEN_RUNNING
               , std::enable_if_t<
                   ! internal::is_Point_set_or_Range_or_Iterator<
                       typename boost::range_value<GraphRange>::type>::value>* = nullptr
#endif
               )
{
  typedef typename boost::range_value<GraphRange>::type                         FaceGraph;
  typedef typename boost::property_map<FaceGraph, boost::vertex_point_t>::type  VPM;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor            vertex_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor              face_descriptor;

  // @todo `Triangle` ought to be just array<int, 3>
  typedef std::vector<int>                                                      Triangle;
  typedef std::vector<Triangle>                                                 TriangleRange;
  typedef std::vector<Point>                                                    PointRange;

  std::vector<PointRange> all_points;
  std::vector<TriangleRange> all_triangles;

  for(const FaceGraph& g : gs)
  {
    PointRange points;
    points.reserve(num_vertices(g));
    TriangleRange triangles;
    triangles.reserve(num_faces(g));

    VPM vpm = get(boost::vertex_point, g);

    // @todo dynamic pmap
    std::unordered_map<typename boost::graph_traits<FaceGraph>::vertex_descriptor, int> vertex_id_map;

    int i = 0;
    for(const vertex_descriptor v : vertices(g))
    {
      points.push_back(get(vpm, v));
      vertex_id_map[v] = i++;
    }

    all_points.push_back(points);
    for(const face_descriptor f : faces(g))
    {
      Triangle triangle;
      for(vertex_descriptor vert : CGAL::vertices_around_face(halfedge(f, g), g))
        triangle.push_back(vertex_id_map[vert]);

      CGAL_assertion(triangle.size() == 3);
      triangles.push_back(triangle);
    }

    all_triangles.push_back(triangles);
  }

  return write_3MF(filename, all_points, all_triangles, names);
}

} } // namespace CGAL::IO

#endif // defined(CGAL_LINKED_WITH_3MF) || defined(DOXYGEN_RUNNING)

#endif // CGAL_BGL_IO_3MF_H
