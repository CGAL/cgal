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

#ifndef CGAL_KSR_3_POLYGON_SPLITTER_H
#define CGAL_KSR_3_POLYGON_SPLITTER_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <queue>

#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Data_structure.h>

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace CGAL
{

namespace KSR_3
{

template <typename GeomTraits>
class Polygon_splitter
{
  typedef typename GeomTraits::Point_2 Point_2;
  typedef typename GeomTraits::Point_3 Point_3;
  
  typedef KSR_3::Data_structure<GeomTraits> Data;
  typedef typename Data::Support_plane Support_plane;
  typedef typename Support_plane::Mesh Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Edge_index Edge_index;
  typedef typename Mesh::Halfedge_index Halfedge_index;
  typedef typename Mesh::Face_index Face_index;
  typedef typename Data::Intersection_graph Intersection_graph;
  typedef typename Intersection_graph::Vertex_descriptor Intersection_vertex;
  typedef typename Intersection_graph::Edge_descriptor Intersection_edge;

  struct Vertex_info
  {
    Vertex_index index;
    Intersection_vertex intersection;

    Vertex_info()
      : index (Vertex_index())
      , intersection (Intersection_graph::null_vertex())
    { }
  };

  struct Face_info
  {
    KSR::size_t index;

    Face_info()
      : index (KSR::uninitialized())
    { }
  };

  typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, GeomTraits> Vbase;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, GeomTraits> Fbase;
  typedef CGAL::Constrained_triangulation_face_base_2<GeomTraits, Fbase> CFbase;
  typedef CGAL::Triangulation_data_structure_2<Vbase, CFbase> TDS;
  typedef CGAL::Constrained_Delaunay_triangulation_2<GeomTraits, TDS, CGAL::Exact_predicates_tag> CDT;
  typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTP;
  typedef typename CDTP::Vertex_handle Vertex_handle;
  typedef typename CDTP::Edge Edge;
  typedef typename CDTP::Face_handle Face_handle;
  typedef typename CDTP::Edge_circulator Edge_circulator;
  typedef typename CDTP::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename CDTP::Finite_faces_iterator Finite_faces_iterator;
  typedef typename CDTP::Constraint_id Cid;
  typedef typename CDTP::Context Context;
  typedef typename CDTP::Context_iterator Context_iterator;
  typedef typename CDTP::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
  
  Data& m_data;
  CDTP m_cdt;
  std::map<Cid, Intersection_edge> m_map_intersections;
  
public:

  Polygon_splitter (Data& data) : m_data (data) { }

  void split_support_plane (KSR::size_t support_plane_idx)
  {
    const Support_plane& support_plane = m_data.support_plane(support_plane_idx);
    Mesh& mesh = m_data.mesh (support_plane_idx);

    // First, insert polygons
    for (Vertex_index vi : mesh.vertices())
    {
      Vertex_handle vh = m_cdt.insert (mesh.point(vi));
      vh->info().index = vi;
    }

    for (Face_index fi : mesh.faces())
    {
      std::vector<Point_2> points;
      for (Halfedge_index hi : halfedges_around_face (halfedge(fi, mesh), mesh))
        points.push_back (mesh.point(mesh.target(hi)));
      points.push_back (points.front());
      Cid cid = m_cdt.insert_constraint (points.begin(), points.end());
      m_map_intersections.insert (std::make_pair (cid, Intersection_graph::null_edge()));
    }

    // Then, add intersection vertices + constraints
    for (const Intersection_edge& intersection_edge : support_plane.intersection_edges())
    {
      Intersection_vertex source = m_data.source(intersection_edge);
      Intersection_vertex target = m_data.target(intersection_edge);
      
      Vertex_handle vsource = m_cdt.insert (support_plane.to_2d(m_data.point_3(source)));
      vsource->info().intersection = source;
      Vertex_handle vtarget = m_cdt.insert (support_plane.to_2d(m_data.point_3(target)));
      vtarget->info().intersection = target;

      Cid cid = m_cdt.insert_constraint (vsource, vtarget);
      m_map_intersections.insert (std::make_pair (cid, intersection_edge));
    }

    // Tag external faces
    std::queue<Face_handle> todo;
    todo.push (m_cdt.incident_faces (m_cdt.infinite_vertex()));
    while (!todo.empty())
    {
      Face_handle fh = todo.front();
      todo.pop();

      if (fh->info().index != KSR::uninitialized())
        continue;

      fh->info().index = KSR::no_element();

      for (int i = 0; i < 3; ++ i)
      {
        Face_handle next = fh->neighbor(i);
        bool border = is_border (std::make_pair(fh, i));

        if (!border)
          todo.push(next);
      }
    }

    KSR::size_t face_index = 0;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
    {
      if (it->info().index != KSR::uninitialized())
        continue;

      todo.push(it);
      KSR::size_t nb_faces = 0;
      while (!todo.empty())
      {
        Face_handle fh = todo.front();
        todo.pop();

        if (fh->info().index != KSR::uninitialized())
          continue;

        fh->info().index = face_index;
        ++ nb_faces;

        for (int i = 0; i < 3; ++ i)
        {
          Face_handle next = fh->neighbor(i);
          bool is_constrained = m_cdt.is_constrained (std::make_pair(fh, i));
          if (!is_constrained)
            todo.push(next);
        }
      }

      ++ face_index;
    }

//    dump(support_plane_idx);

    // Rebuild mesh
    for (Face_index fi : mesh.faces())
      mesh.remove_face(fi);
    for (Edge_index ei : mesh.edges())
      mesh.remove_edge(ei);
    for (Vertex_index vi : mesh.vertices())
      mesh.set_halfedge(vi, Halfedge_index());

    std::set<KSR::size_t> done;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
    {
      CGAL_assertion (it->info().index != KSR::uninitialized());

      if (it->info().index == KSR::no_element())
        continue;

      Edge edge;
      for (int i = 0; i < 3; ++ i)
      {
        edge = std::make_pair (it, i);
        if (m_cdt.is_constrained(edge))
          break;
      }

      if (!m_cdt.is_constrained(edge))
        continue;

      if (!done.insert (edge.first->info().index).second)
        continue;

      std::vector<Vertex_index> new_vertices;

      Edge current = edge;
      do
      {
        Face_handle face = current.first;
        int idx = current.second;
        
        Vertex_handle source = face->vertex (m_cdt.ccw(idx));
        Vertex_handle target = face->vertex (m_cdt.cw(idx));
        if (source->info().index == Vertex_index())
          source->info().index = mesh.add_vertex (source->point());

        new_vertices.push_back (source->info().index);

        Edge next = std::make_pair (face, m_cdt.ccw(idx));
        while (!m_cdt.is_constrained (next))
        {
          Face_handle next_face = next.first->neighbor(next.second);
          CGAL_assertion (next_face->info().index == edge.first->info().index);
          
          int next_idx = m_cdt.ccw(next_face->index (next.first));
          next = std::make_pair (next_face, next_idx);
        }
        CGAL_assertion (next.first->vertex(m_cdt.ccw(next.second)) == target);
        current = next;
      }
      while (current != edge);

      Face_index fi = mesh.add_face (new_vertices);

      CGAL_assertion (fi != Face_index());
    }

    // Set intersection adjacencies
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin();
         it != m_cdt.finite_vertices_end(); ++ it)
      if (it->info().index != Vertex_index()
          && it->info().intersection != Intersection_graph::null_vertex())
      {
        m_data.connect (support_plane_idx, it->info().index, it->info().intersection);
      }

    for (const std::pair<Cid, Intersection_edge>& m : m_map_intersections)
    {
      if (m.second == Intersection_graph::null_edge())
        continue;
      
      Vertices_in_constraint_iterator it = m_cdt.vertices_in_constraint_begin (m.first);
      while (true)
      {
        Vertices_in_constraint_iterator next = it;
        ++ next;
        if (next == m_cdt.vertices_in_constraint_end (m.first))
          break;

        Vertex_handle a = *it;
        Vertex_handle b = *next;

        it = next;

        if (a->info().index == Vertex_index() || b->info().index == Vertex_index())
          continue;

        m_data.connect (support_plane_idx, a->info().index, b->info().index, m.second);
      }
    }
  }

  

private:

  bool is_border (const std::pair<Face_handle, int>& edge) const
  {
    if (!m_cdt.is_constrained(edge))
      return false;

    for (Context_iterator it = m_cdt.contexts_begin (edge.first->vertex((edge.second + 1)%3),
                                                     edge.first->vertex((edge.second + 2)%3));
         it != m_cdt.contexts_end (edge.first->vertex((edge.second + 1)%3),
                                   edge.first->vertex((edge.second + 2)%3)); ++ it)
    {
      typename std::map<Cid, Intersection_edge>::const_iterator
        iter = m_map_intersections.find (it->id());
      if (iter == m_map_intersections.end())
        continue;

      if (iter->second == Intersection_graph::null_edge())
        return true;
    }
    return false;
  }

  void dump(KSR::size_t support_plane_idx)
  {
    typedef CGAL::Surface_mesh<Point_3> Mesh_3;
    typedef typename Mesh_3::template Property_map<typename Mesh_3::Face_index, unsigned char> Uchar_map;

    Mesh_3 mesh;
    Uchar_map red = mesh.template add_property_map<typename Mesh_3::Face_index, unsigned char>("red", 0).first;
    Uchar_map green = mesh.template add_property_map<typename Mesh_3::Face_index, unsigned char>("green", 0).first;
    Uchar_map blue = mesh.template add_property_map<typename Mesh_3::Face_index, unsigned char>("blue", 0).first;

    KSR::size_t bbox_nb_vertices = 0;
    KSR::size_t nb_vertices = 0;

    std::map<Vertex_handle, typename Mesh_3::Vertex_index> map_v2i;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
      map_v2i.insert (std::make_pair
                      (it, mesh.add_vertex (m_data.support_plane(support_plane_idx).to_3d
                                            (it->point()))));

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
    {
      std::array<typename Mesh_3::Vertex_index, 3> vertices;
      for (int i = 0; i < 3; ++ i)
        vertices[i] = map_v2i[it->vertex(i)];
      typename Mesh_3::Face_index face = mesh.add_face (vertices);
      CGAL::Random rand (it->info().index);
      if (it->info().index != KSR::no_element())
      {
        red[face] = (unsigned char)(rand.get_int(32, 192));
        green[face] = (unsigned char)(rand.get_int(32, 192));
        blue[face] = (unsigned char)(rand.get_int(32, 192));
      }
    }
    
    std::string filename = "face_" + std::to_string(support_plane_idx) + ".ply";
    std::ofstream out (filename);
    CGAL::set_binary_mode (out);
    CGAL::write_ply(out, mesh);
  }
  
  
};

}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_POLYGON_SPLITTER_H
