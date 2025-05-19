// Copyright (c) 2019  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Maxime Gimeno

#ifndef CGAL_SURFACE_MESH_IO_3MF_H
#define CGAL_SURFACE_MESH_IO_3MF_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>

#include <CGAL/IO/3MF.h>

#include <iostream>
#include <string>
#include <vector>

#if defined(CGAL_LINKED_WITH_3MF) || defined(DOXYGEN_RUNNING)

#ifdef DOXYGEN_RUNNING
#define CGAL_DEPRECATED
#endif

namespace CGAL {

namespace IO {

// @todo generalize it to any model of `FaceGraph` and put it in BGL/IO (see read_OFF for the face color map)

/*!
 * \ingroup PkgSurfaceMeshIOFunc3MF
 *
 * \brief extracts the surface meshes from an input 3mf file and appends it to `output`.
 *
 * \tparam Point The type of the \em point property of a vertex. There is no requirement on `P`,
 *               besides being default constructible and assignable.
 *               In typical use cases it will be a 2D or 3D point type.
 *
 * \param filename the path to the 3mf file
 * \param output a `std::vector` containing the `CGAL::Surface_mesh`s that will be filled by this function
 *
 * \returns `true` if reading was successful, `false` otherwise.
 */
template<typename Point>
bool read_3MF(const std::string& filename,
              std::vector<CGAL::Surface_mesh<Point> >& output)
{
  typedef std::vector<Point>                                  PointRange;
  typedef std::vector<std::size_t>                            Triangle;
  typedef std::vector<Triangle>                               TriangleRange;
  typedef CGAL::Surface_mesh<Point>                           SMesh;
  typedef typename SMesh::Vertex_index                        Vertex_index;
  typedef typename SMesh::Face_index                          Face_index;

  std::vector<PointRange> all_points;
  std::vector<TriangleRange> all_triangles;
  std::vector<std::string> names;
  std::vector<std::vector<CGAL::IO::Color> > all_colors;

  const bool success = CGAL::IO::read_3MF(filename, all_points, all_triangles, all_colors, names);
  if(!success)
  {
    std::cerr << "Error in reading meshes." << std::endl;
    return false;
  }

  const std::size_t nb_meshes = all_points.size();
  int result = 0;
  output.reserve(nb_meshes);

  for(std::size_t i=0; i<nb_meshes; ++i)
  {
    bool skip = false;
    SMesh sm;
    TriangleRange triangles = all_triangles[i];
    PointRange points = all_points[i];
    std::vector<CGAL::IO::Color> colors = all_colors[i];

    // Create the surface mesh from scratch
    std::size_t n(points.size());
    sm.reserve(n, 0, triangles.size());
    for(const Point& p : points)
      sm.add_vertex(p);

    for(Triangle& triangle : triangles)
    {
      std::vector<Vertex_index> face;
      face.reserve(triangle.size());
      for(std::size_t index : triangle)
        face.push_back(Vertex_index(index));

      Face_index fi = sm.add_face(face);
      if(fi == sm.null_face())
      {
        skip = true;
        sm.clear();
        break;
      }
    }

    if(skip)
      continue;

    const Color& first = colors.front();
    bool need_pmap = false;
    for(const Color& color : colors)
    {
      if(color != first)
      {
        need_pmap = true;
        break;
      }
    }

    if(need_pmap)
    {
      typename SMesh::template Property_map<Face_index, CGAL::IO::Color> fcolor =
          sm.template add_property_map<Face_index,CGAL::IO::Color>("f:color", first).first;

      for(std::size_t pid=0, cs=colors.size(); pid<cs; ++pid)
        put(fcolor, Face_index(pid), colors[pid]); // there can't have any deleted face yet
    }

    output.push_back(sm);
    ++result;
  }

  return true;
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::read_3MF()` should be used instead.
*/
template<typename Point>
CGAL_DEPRECATED int read_3mf(const std::string& filename, std::vector<CGAL::Surface_mesh<Point> >& output)
{
  return IO::read_3MF(filename, output);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // defined(CGAL_LINKED_WITH_3MF) || defined(DOXYGEN_RUNNING)

#endif // CGAL_SURFACE_MESH_IO_3MF_H
