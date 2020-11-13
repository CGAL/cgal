// Copyright (c) 2019 GeometryFactory SARL (France).
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

// #include <CGAL/license/Kinetic_shape_reconstruction.h>

// CGAL includes.
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR_3/Data_structure.h>

namespace CGAL {
namespace KSR_3 {

template<typename GeomTraits>
class Polygon_splitter {

public:
  using Kernel = GeomTraits;

private:
  using FT         = typename Kernel::FT;
  using Point_2    = typename Kernel::Point_2;
  using Point_3    = typename Kernel::Point_3;
  using Line_2     = typename Kernel::Line_2;
  using Vector_2   = typename Kernel::Vector_2;
  using Triangle_2 = typename Kernel::Triangle_2;

  using Data = KSR_3::Data_structure<Kernel>;

  using PVertex = typename Data::PVertex;
  using PFace   = typename Data::PFace;
  using PEdge   = typename Data::PEdge;

  using IVertex = typename Data::IVertex;
  using IEdge   = typename Data::IEdge;

  struct Vertex_info {
    PVertex pvertex;
    IVertex ivertex;
    Vertex_info() :
    pvertex(Data::null_pvertex()),
    ivertex(Data::null_ivertex())
    { }
  };

  struct Face_info {
    KSR::size_t index;
    Face_info() :
    index(KSR::uninitialized())
    { }
  };

  using VBI = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel>;
  using FBI = CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>;
  using CFB = CGAL::Constrained_triangulation_face_base_2<Kernel, FBI>;
  using TDS = CGAL::Triangulation_data_structure_2<VBI, CFB>;
  using TAG = CGAL::Exact_predicates_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, TAG>;
  using TRI = CGAL::Constrained_triangulation_plus_2<CDT>;
  using CID = typename TRI::Constraint_id;

  using Vertex_handle = typename TRI::Vertex_handle;
  using Face_handle   = typename TRI::Face_handle;
  using Edge          = typename TRI::Edge;

  using Mesh_3       = CGAL::Surface_mesh<Point_3>;
  using Vertex_index = typename Mesh_3::Vertex_index;
  using Face_index   = typename Mesh_3::Face_index;
  using Uchar_map    = typename Mesh_3::template Property_map<Face_index, unsigned char>;

  Data& m_data;
  TRI m_cdt;
  std::set<PVertex> m_input;
  std::map<CID, IEdge> m_map_intersections;

public:
  Polygon_splitter(Data& data) :
  m_data(data)
  { }

  void split_support_plane(const KSR::size_t support_plane_idx) {

    // First, insert polygons.
    for (const auto pvertex : m_data.pvertices(support_plane_idx)) {
      const auto vh = m_cdt.insert(m_data.point_2(pvertex));
      vh->info().pvertex = pvertex;
      m_input.insert(pvertex);
    }

    std::vector< std::vector<Point_2> > original_faces;
    std::vector<KSR::size_t> original_input;
    std::vector<Point_2> original_centroids;

    std::vector<Point_2> points;
    std::vector<Triangle_2> triangles;
    for (const PFace pface : m_data.pfaces(support_plane_idx)) {

      points.clear();
      for (const auto pvertex : m_data.pvertices_of_pface(pface)) {
        points.push_back(m_data.point_2(pvertex));
      }

      original_faces.push_back(points);
      original_input.push_back(m_data.input(pface));

      triangles.clear();
      CDT tri;
      for (const auto& point : points) {
        tri.insert(point);
      }
      triangles.reserve(tri.number_of_faces());

      for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
        triangles.push_back(Triangle_2(
            fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point()));
      }
      const auto centroid = CGAL::centroid(triangles.begin(), triangles.end());
      original_centroids.push_back(centroid);

      points.push_back(points.front());
      const auto cid = m_cdt.insert_constraint(points.begin(), points.end());
      m_map_intersections.insert(std::make_pair(cid, Data::null_iedge()));
    }

    // Then, add intersection vertices + constraints.
    for (const auto& iedge : m_data.iedges(support_plane_idx)) {

      const auto source = m_data.source(iedge);
      const auto target = m_data.target(iedge);

      const auto vsource = m_cdt.insert(m_data.to_2d(support_plane_idx, source));
      vsource->info().ivertex = source;
      const auto vtarget = m_cdt.insert(m_data.to_2d(support_plane_idx, target));
      vtarget->info().ivertex = target;

      const auto cid = m_cdt.insert_constraint(vsource, vtarget);
      m_map_intersections.insert(std::make_pair(cid, iedge));
    }

    // Tag external faces.
    std::queue<Face_handle> todo;
    todo.push(m_cdt.incident_faces(m_cdt.infinite_vertex()));
    while (!todo.empty()) {

      const auto fh = todo.front();
      todo.pop();
      if (fh->info().index != KSR::uninitialized()) {
        continue;
      }
      fh->info().index = KSR::no_element();

      for (int i = 0; i < 3; ++i) {
        const auto next = fh->neighbor(i);
        const bool is_on_border = is_border(std::make_pair(fh, i));
        if (!is_on_border) {
          todo.push(next);
        }
      }
    }

    KSR::size_t face_index = 0;
    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      if (fit->info().index != KSR::uninitialized()) {
        continue;
      }

      todo.push(fit);
      KSR::size_t nb_faces = 0;
      while (!todo.empty()) {
        const auto fh = todo.front();
        todo.pop();
        if (fh->info().index != KSR::uninitialized()) {
          continue;
        }
        fh->info().index = face_index;
        ++nb_faces;

        for (int i = 0; i < 3; ++i) {
          const auto next = fh->neighbor(i);
          const bool is_constrained = m_cdt.is_constrained(std::make_pair(fh, i));
          if (!is_constrained) {
            todo.push(next);
          }
        }
      }
      ++face_index;
    }

    // dump(support_plane_idx);
    m_data.clear_polygon_faces(support_plane_idx);

    std::set<KSR::size_t> done;
    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      CGAL_assertion(fit->info().index != KSR::uninitialized());
      if (fit->info().index == KSR::no_element()) {
        continue;
      }

      Edge edge;
      for (int i = 0; i < 3; ++i) {
        edge = std::make_pair(fit, i);
        if (m_cdt.is_constrained(edge)) {
          break;
        }
      }

      if (!m_cdt.is_constrained(edge)) {
        continue;
      }

      if (!done.insert(edge.first->info().index).second) {
        continue;
      }

      std::vector<PVertex> new_vertices;
      auto current = edge;
      do {
        const auto face = current.first;
        const int idx = current.second;

        const auto source = face->vertex(m_cdt.ccw(idx));
        const auto target = face->vertex(m_cdt.cw(idx));
        if (source->info().pvertex == Data::null_pvertex()) {
          source->info().pvertex = m_data.add_pvertex(
            support_plane_idx, source->point());
        }
        new_vertices.push_back(source->info().pvertex);

        auto next = std::make_pair(face, m_cdt.ccw(idx));
        while (!m_cdt.is_constrained(next)) {

          const auto next_face = next.first->neighbor(next.second);
          CGAL_assertion(next_face->info().index == edge.first->info().index);

          const int next_idx = m_cdt.ccw(next_face->index(next.first));
          next = std::make_pair(next_face, next_idx);
        }
        CGAL_assertion(next.first->vertex(m_cdt.ccw(next.second)) == target);
        current = next;

      } while (current != edge);

      const auto pface = m_data.add_pface(new_vertices);
      CGAL_assertion(pface != PFace());

      std::size_t original_idx = 0;
      if (original_faces.size() != 1) {
        // TODO: locate centroid of the face among the different
        // original faces to recover the input index.
        CGAL_assertion_msg(false, "TODO: FINISH THE CASE WITH THE ONE ORIGINAL FACE IN THE POLYGON SPLITTER!");
      }
      m_data.input(pface) = original_input[original_idx];
    }

    // Set intersection adjacencies.
    for (auto vit = m_cdt.finite_vertices_begin(); vit != m_cdt.finite_vertices_end(); ++vit) {
      if (vit->info().pvertex != Data::null_pvertex() &&
          vit->info().ivertex != Data::null_ivertex()) {

        m_data.connect(vit->info().pvertex, vit->info().ivertex);
      }
    }

    for (const auto& m : m_map_intersections) {
      if (m.second == Data::null_iedge()) {
        continue;
      }

      auto vit = m_cdt.vertices_in_constraint_begin(m.first);
      while (true) {
        auto next = vit;
        ++next;
        if (next == m_cdt.vertices_in_constraint_end(m.first)) {
          break;
        }

        const auto a = *vit;
        const auto b = *next;
        vit = next;

        if (a->info().pvertex == Data::null_pvertex() || b->info().pvertex == Data::null_pvertex()) {
          continue;
        }
        m_data.connect(a->info().pvertex, b->info().pvertex, m.second);
      }
    }

    for (const auto pvertex : m_data.pvertices(support_plane_idx)) {
      bool frozen = false;
      auto iedge = Data::null_iedge();
      std::pair<PVertex, PVertex> neighbors(Data::null_pvertex(), Data::null_pvertex());

      for (const auto pedge : m_data.pedges_around_pvertex(pvertex)) {
        if (m_data.has_iedge(pedge)) {
          if (iedge == Data::null_iedge()) {
            iedge = m_data.iedge(pedge);
          } else {
            frozen = true;
            break;
          }
        } else {
          const auto opposite = m_data.opposite(pedge, pvertex);
          if (neighbors.first == Data::null_pvertex()) {
            neighbors.first = opposite;
          } else {
            CGAL_assertion(neighbors.second == Data::null_pvertex());
            neighbors.second = opposite;
          }
        }
      }

      // Several incident intersections = frozen vertex.
      if (frozen) {
        m_data.direction(pvertex) = CGAL::NULL_VECTOR;
        continue;
      }

      // No intersection incident = keep initial direction.
      if (iedge == Data::null_iedge()) {
        continue;
      }
      m_data.connect(pvertex, iedge);
      CGAL_assertion(neighbors.first != Data::null_pvertex() && neighbors.second != Data::null_pvertex());

      bool first_okay = (m_input.find(neighbors.first) != m_input.end());
      auto latest = pvertex;
      auto current = neighbors.first;
      while (!first_okay) {
        const auto pair = m_data.border_prev_and_next(current);
        auto next = pair.first;
        auto ignored = pair.second;

        if (next == latest) {
          std::swap(next, ignored);
        }
        CGAL_assertion(ignored == latest);

        latest = current;
        current = next;
        if (m_input.find(current) != m_input.end()) {
          neighbors.first = current;
          first_okay = true;
        }
      }

      bool second_okay = (m_input.find(neighbors.second) != m_input.end());
      latest = pvertex;
      current = neighbors.second;
      while (!second_okay) {
        const auto pair = m_data.border_prev_and_next(current);
        auto next = pair.first;
        auto ignored = pair.second;

        if (next == latest) {
          std::swap(next, ignored);
        }
        CGAL_assertion(ignored == latest);

        latest = current;
        current = next;
        if (m_input.find(current) != m_input.end()) {
          neighbors.second = current;
          second_okay = true;
        }
      }

      const Line_2 future_line(
        m_data.point_2(neighbors.first, FT(1)), m_data.point_2(neighbors.second, FT(1)));
      const auto intersection_line = m_data.segment_2(support_plane_idx, iedge).supporting_line();
      const Point_2 inter = KSR::intersection<Point_2>(intersection_line, future_line);
      m_data.direction(pvertex) = Vector_2(m_data.point_2(pvertex, FT(0)), inter);
    }
  }

private:

  const bool is_border(const std::pair<Face_handle, int>& edge) const {

    if (!m_cdt.is_constrained(edge)) {
      return false;
    }

    for (auto cit = m_cdt.contexts_begin(
      edge.first->vertex((edge.second + 1) % 3), edge.first->vertex((edge.second + 2) % 3));
      cit != m_cdt.contexts_end(
      edge.first->vertex((edge.second + 1) % 3), edge.first->vertex((edge.second + 2) % 3));
      ++cit) {

      const auto iter = m_map_intersections.find(cit->id());
      if (iter == m_map_intersections.end()) {
        continue;
      }
      if (iter->second == Data::null_iedge()) {
        return true;
      }
    }
    return false;
  }

  void dump(const KSR::size_t support_plane_idx) {

    Mesh_3 mesh;
    Uchar_map red   = mesh.template add_property_map<Face_index, unsigned char>("red", 0).first;
    Uchar_map green = mesh.template add_property_map<Face_index, unsigned char>("green", 0).first;
    Uchar_map blue  = mesh.template add_property_map<Face_index, unsigned char>("blue", 0).first;

    std::map<Vertex_handle, Vertex_index> map_v2i;
    for (auto vit = m_cdt.finite_vertices_begin(); vit != m_cdt.finite_vertices_end(); ++vit) {
      map_v2i.insert(std::make_pair(
        vit, mesh.add_vertex(m_data.support_plane(support_plane_idx).to_3d(vit->point()))));
    }

    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      std::array<Vertex_index, 3> vertices;
      for (std::size_t i = 0; i < 3; ++i) {
        vertices[i] = map_v2i[fit->vertex(i)];
      }

      const auto face = mesh.add_face(vertices);
      CGAL::Random rand(fit->info().index);
      if (fit->info().index != KSR::no_element()) {
        red[face]   = (unsigned char)(rand.get_int(32, 192));
        green[face] = (unsigned char)(rand.get_int(32, 192));
        blue[face]  = (unsigned char)(rand.get_int(32, 192));
      }
    }

    const std::string filename = "face_" + std::to_string(support_plane_idx) + ".ply";
    std::ofstream out(filename);
    out.precision(20);
    CGAL::write_ply(out, mesh);
    out.close();
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_POLYGON_SPLITTER_H
