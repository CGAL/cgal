// Copyright (c) 2019  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Maxime Gimeno


#ifndef CGAL_SURFACE_MESH_IO_3MF_H
#define CGAL_SURFACE_MESH_IO_3MF_H

#include <iostream>
#include <CGAL/IO/read_3mf.h>
#include <CGAL/Surface_mesh.h>

namespace CGAL{
/*!
 * Extracts the surface meshes from an input 3mf file and appends it to `output`.
 *\tparam  Point the Point type of the output meshes.
 * \param file_name the path to the 3mf file.
 * \param output a `std::vector` containing the `CGAL::Surface_mesh`s that will be filled by this function.
 * \return the number of extracted meshes.
 */

template<typename Point>
int read_3mf(const std::string& file_name,
             std::vector<CGAL::Surface_mesh<Point> >& output)
{
  typedef std::vector<Point> PointRange;
  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef typename SMesh::Vertex_index Vertex_index;
  typedef typename SMesh::Face_index Face_index;

  std::vector<PointRange> all_points;
  std::vector<PolygonRange> all_polygons;
  std::vector<std::string> names;
  std::vector<std::vector<CGAL::Color> > all_colors;
  int result = 0;
  int nb_meshes =
      CGAL::read_triangle_soups_from_3mf(file_name,
                                all_points, all_polygons, all_colors, names);
  if(nb_meshes < 0 )
  {
    std::cerr << "Error in reading meshes."<<std::endl;
    return -1;
  }
  output.reserve(nb_meshes);
  for(int i = 0; i< nb_meshes; ++i)
  {
    bool skip = false;
    SMesh sm;
    PolygonRange triangles = all_polygons[i];
    PointRange points = all_points[i];
    std::vector<CGAL::Color> colors = all_colors[i];
    //Create the surface mesh from scratch
    std::size_t n(points.size());
    sm.reserve(n,0, triangles.size());
    for(const Point& p : points)
    {
      sm.add_vertex(p);
    }

    for(Polygon& triangle : triangles)
    {
      std::vector<Vertex_index> face;
      face.reserve(triangle.size());
      for(auto index : triangle)
      {
        face.push_back(Vertex_index(index));
      }
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
    //end constructin the surface mesh from scratch

    CGAL::Color first = colors.front();
    bool need_pmap = false;
    for(auto color : colors)
    {
      if (color != first)
      {
        need_pmap = true;
        break;
      }
    }
    if(need_pmap)
    {
      typename SMesh::template Property_map<Face_index, CGAL::Color> fcolor =
          sm.template add_property_map<Face_index,CGAL::Color>("f:color",first).first;
      for(std::size_t pid = 0; pid < colors.size(); ++pid)
      {
        put(fcolor, Face_index(pid), colors[pid]);//should work bc mesh is just created and shouldn't have any destroyed face.
      }
    }
    output.push_back(sm);
    ++result;
  }
  return result;
}

}//end CGAL
#endif // CGAL_SURFACE_MESH_3MF_H
