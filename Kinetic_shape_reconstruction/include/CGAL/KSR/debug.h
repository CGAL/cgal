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
void dump_intersection_edges (const DS& data, const std::string& tag = std::string())
{
  std::string filename = (tag != std::string() ? tag + "_" : "") + "intersection_edges.polylines.txt";
  std::ofstream out (filename);
  out.precision(18);

  for (const typename DS::Intersection_edge& edge : data.intersection_edges())
  {
//    out << "2 " << data.segment_3 (edge) << std::endl;

    srand (data.source(edge));
    out << "2 " << data.segment_3 (edge).source() + typename DS::Vector_3(0.01 * rand() / double(RAND_MAX),
                                                                          0.01 * rand() / double(RAND_MAX),
                                                                          0.01 * rand() / double(RAND_MAX));
    
    srand (data.target(edge));
    out << " " << data.segment_3 (edge).target() + typename DS::Vector_3(0.01 * rand() / double(RAND_MAX),
                                                                         0.01 * rand() / double(RAND_MAX),
                                                                         0.01 * rand() / double(RAND_MAX))
        << std::endl;
  }
}

template <typename DS>
void dump_constrained_edges (const DS& data, const std::string& tag = std::string())
{
  std::string filename = (tag != std::string() ? tag + "_" : "") + "constrained_edges.polylines.txt";
  std::ofstream out (filename);
  out.precision(18);

  for (KSR::size_t i = 0; i < data.number_of_meshes(); ++ i)
  {
    const typename DS::Mesh& m = data.mesh(i);

    for (const typename DS::Edge_index ei : m.edges())
      if (data.support_plane(i).has_intersection_edge(ei))
        out << "2 " << data.support_plane(i).segment_3 (ei, 0) << std::endl;
  }
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

  KSR::size_t bbox_nb_vertices = 0;
  KSR::size_t nb_vertices = 0;

  KSR::vector<typename Mesh::Vertex_index> vertices;
  for (KSR::size_t i = 0; i < data.number_of_meshes(); ++ i)
  {
    const typename DS::Mesh& m = data.mesh(i);

    if (data.is_bbox_mesh(i))
    {
      KSR::size_t new_vertices = 0;
      for (typename DS::Vertex_index vi : m.vertices())
      {
        bbox_mesh.add_vertex (data.point_of_vertex (i, vi));
        ++ new_vertices;
      }
      
      for (typename DS::Face_index fi : m.faces())
      {
        vertices.clear();
        for(typename DS::Halfedge_index hi : halfedges_around_face(halfedge(fi, m),m))
          vertices.push_back (typename Mesh::Vertex_index(KSR::size_t(source(hi, m)) + bbox_nb_vertices));
        
        typename Mesh::Face_index face = bbox_mesh.add_face (vertices);
        std::tie (bbox_red[face], bbox_green[face], bbox_blue[face])
          = get_idx_color ((i+1) * (fi+1));
      }
      bbox_nb_vertices += new_vertices;
    }
    else
    {
      KSR::size_t new_vertices = 0;
      for (typename DS::Vertex_index vi : m.vertices())
      {
        mesh.add_vertex (data.point_of_vertex (i, vi));
        ++ new_vertices;
      }
      
      for (typename DS::Face_index fi : m.faces())
      {
        vertices.clear();
        for(typename DS::Halfedge_index hi : halfedges_around_face(halfedge(fi, m),m))
          vertices.push_back (typename Mesh::Vertex_index(KSR::size_t(source(hi, m)) + nb_vertices));
        
        typename Mesh::Face_index face = mesh.add_face (vertices);
        std::tie (red[face], green[face], blue[face])
          = get_idx_color (i * (fi+1));
      }
      nb_vertices += new_vertices;
    }
  }
    
  std::string filename = (tag != std::string() ? tag + "_" : "") + "polygons.ply";
  std::ofstream out (filename);
  CGAL::set_binary_mode (out);
  CGAL::write_ply(out, mesh);

#if 1
  std::string bbox_filename = (tag != std::string() ? tag + "_" : "") + "bbox_polygons.ply";
  std::ofstream bbox_out (bbox_filename);
  CGAL::set_binary_mode (bbox_out);
  CGAL::write_ply(bbox_out, bbox_mesh);
#endif
  
}

template <typename DS, typename Event>
void dump_event (const DS& data, const Event& ev, const std::string& tag = std::string())
{
  std::string lfilename = (tag != std::string() ? tag + "_" : "") + "event_line.polylines.txt";
  std::ofstream lout (lfilename);
  lout.precision(18);

  // TODO

  std::string vfilename = (tag != std::string() ? tag + "_" : "") + "event_vertex.xyz";
  std::ofstream vout (vfilename);
  vout.precision(18);
  vout << data.point_of_vertex(ev.vertex_idx()) << std::endl;
}

template <typename DS>
void dump (const DS& data, const std::string& tag = std::string())
{
  dump_intersection_edges (data, tag);
  dump_constrained_edges (data, tag);
  dump_polygons (data, tag);
}


} // namespace KSR_3
} // namespace CGAL



#endif // CGAL_KSR_DEBUG_H
