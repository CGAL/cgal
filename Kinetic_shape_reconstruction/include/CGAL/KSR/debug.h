// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_DEBUG_H
#define CGAL_KSR_DEBUG_H

#include <CGAL/KSR/utils.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Random.h>

namespace CGAL
{
namespace KSR_3
{

std::tuple<unsigned char, unsigned char, unsigned char>
get_idx_color (KSR::size_t idx)
{
  CGAL::Random rand (idx);
  return std::make_tuple ((unsigned char)(rand.get_int(32, 192)),
                          (unsigned char)(rand.get_int(32, 192)),
                          (unsigned char)(rand.get_int(32, 192)));
}

template <typename DS>
void dump_vertices (const DS& data, const std::string& tag = std::string())
{
  typedef CGAL::Point_set_3<typename DS::Point_3> Point_set;
  typedef typename Point_set::template Property_map<unsigned char> Uchar_map;
  
  Point_set points;
  points.add_normal_map();
  Uchar_map red = points.template add_property_map<unsigned char>("red", 0).first;
  Uchar_map green = points.template add_property_map<unsigned char>("green", 0).first;
  Uchar_map blue = points.template add_property_map<unsigned char>("blue", 0).first;

  for (KSR::size_t i = 0; i < data.number_of_vertices(); ++ i)
  {
    typename Point_set::iterator it
      = points.insert (data.point_of_vertex(i), data.direction_of_vertex(i));
    std::tie (red[*it], green[*it], blue[*it])
      = get_idx_color (data.vertex(i).polygon_idx());
  }

  std::string filename = (tag != std::string() ? tag + "_" : "") + "vertices.ply";
  std::ofstream out (filename);
  CGAL::set_binary_mode (out);
  out << points;
}

template <typename DS>
void dump_intersection_lines (const DS& data, const std::string& tag = std::string())
{
  std::string filename = (tag != std::string() ? tag + "_" : "") + "intersection_lines.polylines.txt";
  std::ofstream out (filename);
  out.precision(18);
  
  for (KSR::size_t i = 0; i < data.number_of_intersection_lines(); ++ i)
    out << "2 "
        << data.meta_vertex (data.intersection_line(i).meta_vertices_idx()[0]).point() << " "
        << data.meta_vertex (data.intersection_line(i).meta_vertices_idx()[1]).point() << std::endl;
}

template <typename DS>
void dump_segments (const DS& data, const std::string& tag = std::string())
{
  std::string filename = (tag != std::string() ? tag + "_" : "") + "segments.polylines.txt";
  std::ofstream out (filename);
  out.precision(18);
  
  for (KSR::size_t i = 0; i < data.number_of_segments(); ++ i)
    out << "2 "
        << data.point_of_vertex (data.segment(i).source_idx()) << " "
        << data.point_of_vertex (data.segment(i).target_idx()) << std::endl;
}

template <typename DS>
void dump_polygons (const DS& data, const std::string& tag = std::string())
{
  typedef CGAL::Surface_mesh<typename DS::Point_3> Mesh;
  typedef typename Mesh::template Property_map<typename Mesh::Face_index, unsigned char> Uchar_map;

  Mesh mesh;
  Uchar_map red = mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("red", 0).first;
  Uchar_map green = mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("green", 0).first;
  Uchar_map blue = mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("blue", 0).first;
  
  Mesh bbox_mesh;
  Uchar_map bbox_red = bbox_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("red", 0).first;
  Uchar_map bbox_green = bbox_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("green", 0).first;
  Uchar_map bbox_blue = bbox_mesh.template add_property_map<typename Mesh::Face_index, unsigned char>("blue", 0).first;
  
  std::vector<typename Mesh::Vertex_index> vertices;
  for (KSR::size_t i = 0; i < data.number_of_polygons(); ++ i)
  {
    vertices.clear();

    if (data.is_bbox_polygon(i))
    {
      for (KSR::size_t vertex_idx : data.polygon(i).vertices_idx())
        vertices.push_back (bbox_mesh.add_vertex (data.point_of_vertex(vertex_idx)));
      typename Mesh::Face_index idx = bbox_mesh.add_face(vertices);
      std::tie (bbox_red[idx], bbox_green[idx], bbox_blue[idx])
      = get_idx_color (i);
    }
    else
    {
      for (KSR::size_t vertex_idx : data.polygon(i).vertices_idx())
        vertices.push_back (mesh.add_vertex (data.point_of_vertex(vertex_idx)));
      typename Mesh::Face_index idx = mesh.add_face(vertices);
      std::tie (red[idx], green[idx], blue[idx])
        = get_idx_color (i);
    }
  }

  std::string filename = (tag != std::string() ? tag + "_" : "") + "polygons.ply";
  std::ofstream out (filename);
  CGAL::set_binary_mode (out);
  CGAL::write_ply(out, mesh);
  
  std::string bbox_filename = (tag != std::string() ? tag + "_" : "") + "bbox_polygons.ply";
  std::ofstream bbox_out (bbox_filename);
  CGAL::set_binary_mode (bbox_out);
  CGAL::write_ply(bbox_out, bbox_mesh);

}

template <typename DS>
void dump_meta_vertices (const DS& data, const std::string& tag = std::string())
{
  typedef CGAL::Point_set_3<typename DS::Point_3> Point_set;
  
  Point_set points;

  for (KSR::size_t i = 0; i < data.number_of_meta_vertices(); ++ i)
    points.insert (data.meta_vertex(i).point());

  std::string filename = (tag != std::string() ? tag + "_" : "") + "meta_vertices.ply";
  std::ofstream out (filename);
  CGAL::set_binary_mode (out);
  out << points;
}

template <typename DS>
void dump (const DS& data, const std::string& tag = std::string())
{
  dump_vertices (data, tag);
  dump_intersection_lines (data, tag);
  dump_segments (data, tag);
  dump_polygons (data, tag);
  dump_meta_vertices (data, tag);
}


} // namespace KSR_3
} // namespace CGAL



#endif // CGAL_KSR_DEBUG_H
