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
  typedef typename GeomTraits::Segment_2 Segment_2;
  typedef typename GeomTraits::Line_2 Line_2;
  typedef typename GeomTraits::Vector_2 Vector_2;
  
  typedef KSR_3::Data_structure<GeomTraits> Data;
  typedef typename Data::PVertex PVertex;
  typedef typename Data::PEdge PEdge;
  typedef typename Data::PFace PFace;
  typedef typename Data::IEdge IEdge;
  typedef typename Data::IVertex IVertex;

  struct Vertex_info
  {
    PVertex pvertex;
    IVertex ivertex;

    Vertex_info()
      : pvertex (Data::null_pvertex())
      , ivertex (Data::null_ivertex())
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
  std::map<Cid, IEdge> m_map_intersections;
  
public:

  Polygon_splitter (Data& data) : m_data (data) { }

  void split_support_plane (KSR::size_t support_plane_idx, unsigned int k)
  {
    // First, insert polygons
    for (PVertex pvertex : m_data.pvertices(support_plane_idx))
    {
      Vertex_handle vh = m_cdt.insert (m_data.point_2 (pvertex));
      vh->info().pvertex = pvertex;
    }

    std::vector<std::vector<Point_2> > original_faces;
    std::vector<KSR::size_t> original_input;
    std::vector<Point_2> original_centroids;

    for (PFace pface : m_data.pfaces(support_plane_idx))
    {
      std::vector<Point_2> points;
      for (PVertex pvertex : m_data.pvertices_of_pface (pface))
        points.push_back (m_data.point_2(pvertex));
      
      original_faces.push_back (points);
      original_input.push_back (m_data.input(pface));
      original_centroids.push_back
        (CGAL::centroid (points.begin(), points.end()));
      
      points.push_back (points.front());
      Cid cid = m_cdt.insert_constraint (points.begin(), points.end());
      m_map_intersections.insert (std::make_pair (cid, Data::null_iedge()));
    }

    // Then, add intersection vertices + constraints
    for (const IEdge& iedge : m_data.iedges(support_plane_idx))
    {
      IVertex source = m_data.source(iedge);
      IVertex target = m_data.target(iedge);
      
      Vertex_handle vsource = m_cdt.insert (m_data.to_2d(support_plane_idx, source));
      vsource->info().ivertex = source;
      Vertex_handle vtarget = m_cdt.insert (m_data.to_2d(support_plane_idx, target));
      vtarget->info().ivertex = target;

      Cid cid = m_cdt.insert_constraint (vsource, vtarget);
      m_map_intersections.insert (std::make_pair (cid, iedge));
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

    m_data.clear_polygon_faces (support_plane_idx);

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

      std::vector<PVertex> new_vertices;

      Edge current = edge;
      do
      {
        Face_handle face = current.first;
        int idx = current.second;
        
        Vertex_handle source = face->vertex (m_cdt.ccw(idx));
        Vertex_handle target = face->vertex (m_cdt.cw(idx));
        if (source->info().pvertex == Data::null_pvertex())
          source->info().pvertex = m_data.add_pvertex (support_plane_idx, source->point());

        new_vertices.push_back (source->info().pvertex);

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

      PFace pface = m_data.add_pface (new_vertices);
      CGAL_assertion (pface != PFace());

      m_data.k(pface) = k;

      std::size_t original_idx = 0;
      if (original_faces.size() != 1)
      {
        // TODO: locate centroid of the face among the different
        // original faces to recover the input index
        CGAL_assertion_msg(false, "TODO!");
      }

      m_data.input(pface) = original_input[original_idx];
      for (PVertex pvertex : new_vertices)
        m_data.direction(pvertex) = KSR::normalize (Vector_2 (original_centroids[original_idx],
                                                              m_data.point_2(pvertex)));
    }

    // Set intersection adjacencies
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin();
         it != m_cdt.finite_vertices_end(); ++ it)
      if (it->info().pvertex != Data::null_pvertex()
          && it->info().ivertex != Data::null_ivertex())
      {
        m_data.connect (it->info().pvertex, it->info().ivertex);
      }

    for (const std::pair<Cid, IEdge>& m : m_map_intersections)
    {
      if (m.second == Data::null_iedge())
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

        if (a->info().pvertex == Data::null_pvertex() || b->info().pvertex == Data::null_pvertex())
          continue;

        m_data.connect (a->info().pvertex, b->info().pvertex, m.second);
      }
    }


    for (const PVertex pvertex : m_data.pvertices(support_plane_idx))
    {
      bool frozen = false;
      IEdge iedge = Data::null_iedge();

      std::pair<PVertex, PVertex> neighbors (Data::null_pvertex(), Data::null_pvertex());
      
      for (PEdge pedge : m_data.pedges_around_pvertex (pvertex))
      {
        if (m_data.has_iedge (pedge))
        {
          if (iedge == Data::null_iedge())
            iedge = m_data.iedge(pedge);
          else
          {
            frozen = true;
            break;
          }
        }
        else
        {
          PVertex opposite = m_data.opposite (pedge, pvertex);
          if (neighbors.first == Data::null_pvertex())
            neighbors.first = opposite;
          else
          {
            CGAL_assertion (neighbors.second == Data::null_pvertex());
            neighbors.second = opposite;
          }
        }
      }

      // Several incident intersections = frozen vertex
      if (frozen)
      {
        m_data.direction(pvertex) = CGAL::NULL_VECTOR;
        continue;
      }

      // No intersection incident = keep initial direction
      if (iedge == Data::null_iedge())
        continue;

      m_data.connect (pvertex, iedge);
      
      CGAL_assertion (neighbors.first != Data::null_pvertex() && neighbors.second != Data::null_pvertex());

      Line_2 future_line (m_data.point_2 (neighbors.first, 1),
                          m_data.point_2 (neighbors.second, 1));

      Line_2 intersection_line = m_data.segment_2 (support_plane_idx, iedge).supporting_line();

      Point_2 inter = KSR::intersection_2<Point_2> (intersection_line, future_line);

      m_data.direction(pvertex) = Vector_2 (m_data.point_2(pvertex, 0), inter);
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
      typename std::map<Cid, IEdge>::const_iterator
        iter = m_map_intersections.find (it->id());
      if (iter == m_map_intersections.end())
        continue;

      if (iter->second == Data::null_iedge())
        return true;
    }
    return false;
  }

  // void dump(KSR::size_t support_plane_idx)
  // {
  //   typedef CGAL::Surface_mesh<Point_3> Mesh_3;
  //   typedef typename Mesh_3::template Property_map<typename Mesh_3::Face_index, unsigned char> Uchar_map;

  //   Mesh_3 mesh;
  //   Uchar_map red = mesh.template add_property_map<typename Mesh_3::Face_index, unsigned char>("red", 0).first;
  //   Uchar_map green = mesh.template add_property_map<typename Mesh_3::Face_index, unsigned char>("green", 0).first;
  //   Uchar_map blue = mesh.template add_property_map<typename Mesh_3::Face_index, unsigned char>("blue", 0).first;

  //   KSR::size_t bbox_nb_vertices = 0;
  //   KSR::size_t nb_vertices = 0;

  //   std::map<Vertex_handle, typename Mesh_3::Vertex_index> map_v2i;
  //   for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
  //     map_v2i.insert (std::make_pair
  //                     (it, mesh.add_vertex (m_data.support_plane(support_plane_idx).to_3d
  //                                           (it->point()))));

  //   for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
  //   {
  //     std::array<typename Mesh_3::Vertex_index, 3> vertices;
  //     for (int i = 0; i < 3; ++ i)
  //       vertices[i] = map_v2i[it->vertex(i)];
  //     typename Mesh_3::Face_index face = mesh.add_face (vertices);
  //     CGAL::Random rand (it->info().index);
  //     if (it->info().index != KSR::no_element())
  //     {
  //       red[face] = (unsigned char)(rand.get_int(32, 192));
  //       green[face] = (unsigned char)(rand.get_int(32, 192));
  //       blue[face] = (unsigned char)(rand.get_int(32, 192));
  //     }
  //   }
    
  //   std::string filename = "face_" + std::to_string(support_plane_idx) + ".ply";
  //   std::ofstream out (filename);
  //   CGAL::set_binary_mode (out);
  //   CGAL::write_ply(out, mesh);
  // }
  
  
};

}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_POLYGON_SPLITTER_H
