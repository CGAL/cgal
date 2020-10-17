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
  // Kernel types.
  using Kernel    = GeomTraits;
  using FT        = typename Kernel::FT;
  using Point_2   = typename Kernel::Point_2;
  using Point_3   = typename Kernel::Point_3;
  using Vector_2  = typename Kernel::Vector_2;
  using Segment_2 = typename Kernel::Segment_2;
  using Line_2    = typename Kernel::Line_2;

  // Data structure types.
  using Data_structure = KSR_3::Data_structure<Kernel>;
  using Support_plane  = typename Data_structure::Support_plane;
  using IGraph         = typename Data_structure::IGraph;

  // Support plane types.
  using Vertex_index = typename Support_plane::Vertex_index;
  using PVertex      = typename Support_plane::PVertex;
  using PEdge        = typename Support_plane::PEdge;
  using PFace        = typename Support_plane::PFace;

  // Intersection graph types.
  using IEdge   = typename IGraph::Edge_descriptor;
  using IVertex = typename IGraph::Vertex_descriptor;

  using Saver = KSR_3::Saver<Kernel>;

  // Triangulation vertex info.
  struct Vertex_info {
    PVertex pvertex;
    IVertex ivertex;
    Vertex_info() :
      pvertex(Support_plane::null_pvertex()),
      ivertex(IGraph::null_vertex())
    { }
  };

  // Triangulation face info.
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

public:

  Polygon_splitter(Data_structure& data) :
    m_data(data)
  { }

  void split_support_plane(
    const KSR::size_t support_plane_idx, const unsigned int k) {

    auto& support_plane = m_data.support_planes()[support_plane_idx];
    KSR::vector< KSR::vector<Point_2> > original_faces;
    KSR::vector<KSR::size_t> original_input;
    KSR::vector<Point_2> original_centroids;

    // Create CDT.
    create_triangulation(
      support_plane_idx,
      original_faces, original_input, original_centroids);
    add_igraph_edges(
      support_plane_idx);
    tag_faces();

    // dump_cdt("debug-tagged-cdt-sp-" + std::to_string(support_plane_idx));

    // Remove all original faces from the support plane.
    support_plane.clear_faces();

    // Split polygons.
    create_new_faces(
      original_faces, original_input, original_centroids,
      support_plane_idx, k);
    set_intersection_adjacencies(support_plane_idx);
    set_intersection_directions(support_plane_idx);
  }

private:
  Data_structure& m_data;
  TRI m_cdt;
  KSR::set<Vertex_index> m_input;
  KSR::map<CID, IEdge> m_map_intersections;
  const Saver m_saver;

  // Create triangulation based on all input polygons from the current support plane.
  void create_triangulation(
    const KSR::size_t support_plane_idx,
    KSR::vector< KSR::vector<Point_2> >& original_faces,
    KSR::vector<KSR::size_t>& original_input,
    KSR::vector<Point_2>& original_centroids) {

    m_cdt.clear();
    original_faces.clear();
    original_input.clear();
    original_centroids.clear();

    const auto current_time = m_data.current_time();
    CGAL_assertion(current_time == FT(0));

    const auto& support_plane = m_data.support_planes()[support_plane_idx];

    // Insert CDT vertices.
    for (const auto pvertex : support_plane.pvertices()) {
      const auto vertex_index = pvertex.second;
      const auto vertex = support_plane.point_2(vertex_index, current_time);
      const auto vh = m_cdt.insert(vertex);
      vh->info().pvertex = pvertex;
      m_input.insert(vertex_index);
    }

    // Insert CDT constraints for all input polygons.
    KSR::vector<Point_2> vertices;
    for (const auto pface : support_plane.pfaces()) {
      const auto face_index = pface.second;

      // Get vertices of the original input polygon.
      vertices.clear();
      for (const auto pvertex : support_plane.pvertices_of_pface(pface)) {
        const auto vertex_index = pvertex.second;
        const auto vertex = support_plane.point_2(vertex_index, current_time);
        vertices.push_back(vertex);
      }

      // Add face and the corresponding index of the input polygon.
      original_faces.push_back(vertices);
      original_input.push_back(support_plane.input(face_index));
      original_centroids.push_back(
        CGAL::centroid(vertices.begin(), vertices.end()));

      // Insert constraints. TODO: Should we use vertex handles here instead?
      vertices.push_back(vertices.front()); // close this face by repeating the first vertex

      // TODO: Is this cid for the whole polygon boundary?
      const CID cid = m_cdt.insert_constraint(vertices.begin(), vertices.end());
      m_map_intersections.insert(std::make_pair(cid, IGraph::null_edge()));
    }
  }

  // Add intersection vertices and constraints for all iedges, which belong
  // to the current support plane.
  void add_igraph_edges(const KSR::size_t support_plane_idx) {

    const auto& igraph = m_data.igraph();
    const auto& support_plane = m_data.support_planes()[support_plane_idx];
    const auto& iedges = support_plane.data().iedges;

    for (const auto& iedge : iedges) {
      const auto source = igraph.source(iedge);
      const auto target = igraph.target(iedge);

      const auto vh_source = m_cdt.insert(m_data.point_2(support_plane_idx, source));
      vh_source->info().ivertex = source;
      const auto vh_target = m_cdt.insert(m_data.point_2(support_plane_idx, target));
      vh_target->info().ivertex = target;

      const CID cid = m_cdt.insert_constraint(vh_source, vh_target);
      m_map_intersections.insert(std::make_pair(cid, iedge));
    }
  }

  // Tag all faces in CDT by splitting them into external and internal.
  void tag_faces() {

    KSR::queue<Face_handle> todo;
    todo.push(m_cdt.incident_faces(m_cdt.infinite_vertex())); // start from the infinite face

    // Initialize all external faces by setting KSR::no_element().
    while (!todo.empty()) {
      const auto fh = todo.front();
      todo.pop();

      if (fh->info().index != KSR::uninitialized()) // skip those, which are already handled
        continue;
      fh->info().index = KSR::no_element(); // setting face index

      for (KSR::size_t i = 0; i < 3; ++i) {
        const auto next = fh->neighbor(i);
        const auto edge = std::make_pair(fh, i);
        const bool is_border_edge = is_border(edge); // border between external and internal faces
        if (!is_border_edge) todo.push(next);
      }
    }

    // Setting indices of all internal faces.
    KSR::size_t face_index = 0;
    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      if (fit->info().index != KSR::uninitialized()) // skip external faces
        continue;

      todo.push(fit);
      while (!todo.empty()) {
        const auto fh = todo.front();
        todo.pop();

        if (fh->info().index != KSR::uninitialized()) // skip those, which are already handled
          continue;
        fh->info().index = face_index; // setting face index

        for (KSR::size_t i = 0; i < 3; ++i) {
          const auto next = fh->neighbor(i);
          const auto edge = std::make_pair(fh, i);
          const bool is_constrained_edge = m_cdt.is_constrained(edge);
          if (!is_constrained_edge)
            todo.push(next);
        }
      }
      ++face_index;
    }
  }

  // Check if this triangulation edge is a border edge between external and internal faces.
  const bool is_border(const std::pair<Face_handle, int>& edge) const {

    if (!m_cdt.is_constrained(edge)) // skip unconstrained edges, they can't be borders
      return false;

    // If it is constrained, then:
    const auto fh = edge.first;
    const auto i  = edge.second;
    const auto ip = (i + 1) % 3; // next neighbor
    const auto im = (i + 2) % 3; // prev neighbor

    for (auto cit = m_cdt.contexts_begin(fh->vertex(ip), fh->vertex(im)); // context of this edge
      cit != m_cdt.contexts_end(fh->vertex(ip), fh->vertex(im)); ++cit) {

      const auto iter = m_map_intersections.find(cit->id()); // find original constraint
      if (iter == m_map_intersections.end()) // if no, skip
        continue;
      if (iter->second == IGraph::null_edge())
        return true; // all input edges are marked as null
    }
    return false;
  }

  // Create new splitted faces and set basic information like vertex directions, k, etc.
  void create_new_faces(
    const KSR::vector< KSR::vector<Point_2> >& original_faces,
    const KSR::vector<KSR::size_t>& original_input,
    const KSR::vector<Point_2>& original_centroids,
    const KSR::size_t support_plane_idx,
    const unsigned int k) {

    const auto current_time = m_data.current_time();
    CGAL_assertion(current_time == FT(0));

    auto& support_plane = m_data.support_planes()[support_plane_idx];
    auto& mesh = support_plane.data().mesh;

    KSR::set<KSR::size_t> done;
    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      CGAL_assertion(fit->info().index != KSR::uninitialized());
      if (fit->info().index == KSR::no_element()) // skip external faces
        continue;

      // Try to find a constrained edge.
      std::pair<Face_handle, KSR::size_t> edge;
      for (KSR::size_t i = 0; i < 3; ++i) {
        edge = std::make_pair(fit, i);
        if (m_cdt.is_constrained(edge))
          break;
      }

      // If no constrained edge found, skip this face.
      if (!m_cdt.is_constrained(edge))
        continue;
      if (!done.insert(edge.first->info().index).second) // if not inserted, skip
        continue;

      // Traverse the boundary of the polygon and create mesh vertices.
      KSR::vector<Vertex_index> new_vertices;
      auto current = edge;
      do {
        const auto face = current.first;
        const auto idx  = current.second;

        // Add first edge.
        const auto source = face->vertex(m_cdt.ccw(idx));
        const auto target = face->vertex(m_cdt.cw(idx));
        if (source->info().pvertex == Support_plane::null_pvertex()) {
          const auto vertex_index = mesh.add_vertex(source->point());
          source->info().pvertex = PVertex(support_plane_idx, vertex_index);
        }
        new_vertices.push_back(source->info().pvertex.second);

        // Find next edge.
        auto next = std::make_pair(face, m_cdt.ccw(idx));
        while (!m_cdt.is_constrained(next)) {

          const auto next_face = next.first->neighbor(next.second);
          CGAL_assertion(next_face->info().index == edge.first->info().index);

          const auto next_idx = m_cdt.ccw(next_face->index(next.first));
          next = std::make_pair(next_face, next_idx);
        }
        CGAL_assertion(next.first->vertex(m_cdt.ccw(next.second)) == target);
        current = next;

      } while (current != edge); // until we come back

      // Add new mesh face.
      const auto face_index = mesh.add_face(new_vertices);
      const PFace pface(support_plane_idx, face_index);
      CGAL_assertion(pface != PFace());

      // Set face number of intersections.
      support_plane.set_k(face_index, k);

      // Set face directions and index of the original face.
      KSR::size_t original_idx = 0;
      if (original_faces.size() != 1) {
        // TODO: locate centroid of the face among the different
        // original faces to recover the input index.
        CGAL_assertion_msg(false, "TODO: FINISH THE CASE WITH THE ONE ORIGINAL FACE IN THE POLYGON SPLITTER!");
      }

      support_plane.set_finput_map(face_index, original_input[original_idx]);
      for (const auto& vertex_index : new_vertices) {
        const auto& centroid = original_centroids[original_idx];
        const auto point = support_plane.point_2(vertex_index, current_time);
        Vector_2 direction(centroid, point);
        KSR::normalize(direction);
        support_plane.set_direction_map(vertex_index, direction);
      }
    }
    // std::cout << "number of new faces: " << mesh.number_of_faces() << std::endl;
  }

  // Set ivertices and iedges of the new face.
  void set_intersection_adjacencies(
    const KSR::size_t support_plane_idx) {

    auto& support_plane = m_data.support_planes()[support_plane_idx];

    // Set ivertices.
    for (auto vit = m_cdt.finite_vertices_begin(); vit != m_cdt.finite_vertices_end(); ++vit) {
      if (vit->info().pvertex != Support_plane::null_pvertex() && vit->info().ivertex != IGraph::null_vertex()) {
        const auto vertex_index = vit->info().pvertex.second;
        support_plane.set_ivertex_map(vertex_index, vit->info().ivertex);
      }
    }

    // Set iedges.
    for (const auto& pair : m_map_intersections) {
      const auto& cid = pair.first;
      const auto& iedge = pair.second;

      if (iedge == IGraph::null_edge())
        continue;

      auto vit = m_cdt.vertices_in_constraint_begin(cid);
      while (true) {

        auto next = vit;
        ++next;
        if (next == m_cdt.vertices_in_constraint_end(cid))
          break;

        auto vh_a = *vit;
        auto vh_b = *next;

        vit = next;
        if (vh_a->info().pvertex == Support_plane::null_pvertex() || vh_b->info().pvertex == Support_plane::null_pvertex())
          continue;

        const auto va = vh_a->info().pvertex.second;
        const auto vb = vh_b->info().pvertex.second;
        support_plane.set_eiedge_map(va, vb, iedge);
      }
    }
  }

  // Set directions of the intersection points or froze them.
  void set_intersection_directions(
    const KSR::size_t support_plane_idx) {

    const auto current_time = m_data.current_time();
    CGAL_assertion(current_time == FT(0));

    auto& support_plane = m_data.support_planes()[support_plane_idx];

    // Go through all mesh vertices.
    for (const auto pvertex : support_plane.pvertices()) {
      const auto vertex_index = pvertex.second;

      bool frozen = false;
      auto iedge = IGraph::null_edge();
      std::pair<Vertex_index, Vertex_index> neighbors;
      neighbors.first  = Vertex_index();
      neighbors.second = Vertex_index();

      // Find neighbors of the vertex.
      for (const auto pedge : support_plane.pedges_around_pvertex(pvertex)) {
        const auto edge_index = pedge.second;

        if (support_plane.has_iedge(edge_index)) {
          if (iedge == IGraph::null_edge()) {
            iedge = support_plane.iedge(edge_index);
          } else {
            frozen = true; // we found a frozen vertex
            break;
          }
        } else {
          const auto opposite = support_plane.opposite(edge_index, vertex_index);
          if (neighbors.first == Vertex_index()) {
            neighbors.first = opposite;
          } else {
            CGAL_assertion(neighbors.second == Vertex_index());
            neighbors.second = opposite;
          }
        }
      }

      // Several incident intersections = frozen vertex.
      if (frozen) {
        support_plane.set_direction_map(vertex_index, CGAL::NULL_VECTOR);
        continue;
      }

      // No intersection incident = keep initial direction.
      if (iedge == IGraph::null_edge())
        continue;

      // Otherwise.
      support_plane.set_viedge_map(vertex_index, iedge);

      // TODO: Why this assertion fails for bbox support planes?
      CGAL_assertion(
        neighbors.first != Vertex_index() && neighbors.second != Vertex_index());

      // Find the first neighbor along the border.
      bool first_okay = (m_input.find(neighbors.first) != m_input.end());
      Vertex_index latest  = vertex_index;
      Vertex_index current = neighbors.first;
      while (!first_okay) {
        const auto pair = support_plane.border_prev_and_next(current);
        auto next    = pair.first;
        auto ignored = pair.second;

        if (next == latest)
          std::swap(next, ignored);
        CGAL_assertion(ignored == latest);

        latest  = current;
        current = next;
        if (m_input.find(current) != m_input.end()) {
          neighbors.first = current;
          first_okay = true;
        }
      }

      // Find the second neighbor along the border.
      bool second_okay = (m_input.find(neighbors.second) != m_input.end());
      latest  = vertex_index;
      current = neighbors.second;
      while (!second_okay) {
        const auto pair = support_plane.border_prev_and_next(current);
        auto next    = pair.first;
        auto ignored = pair.second;

        if (next == latest)
          std::swap(next, ignored);
        CGAL_assertion(ignored == latest);

        latest  = current;
        current = next;
        if (m_input.find(current) != m_input.end()) {
          neighbors.second = current;
          second_okay = true;
        }
      }

      // Set direction.
      const FT future_time = FT(1);
      const Line_2 future_line(
        support_plane.point_2(neighbors.first, future_time),
        support_plane.point_2(neighbors.second, future_time));

      const auto segment = m_data.segment_2(support_plane_idx, iedge);
      const Line_2 intersection_line = segment.supporting_line();
      const Point_2 intersection_point = KSR::intersection<Point_2>(
        intersection_line, future_line);

      // Do not normalize here!
      const auto point = support_plane.point_2(vertex_index, current_time);
      const Vector_2 direction(point, intersection_point);
      support_plane.set_direction_map(vertex_index, direction);
    }
  }

  void dump_cdt(const std::string file_name) const {

    KSR::vector<Point_3> polygon;
    KSR::vector< KSR::vector<Point_3> > polygons;
    KSR::vector<Color> colors;

    const KSR::size_t num_faces = m_cdt.number_of_faces();
    polygons.reserve(num_faces);
    colors.reserve(num_faces);

    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      const auto& p0 = fit->vertex(0)->point();
      const auto& p1 = fit->vertex(1)->point();
      const auto& p2 = fit->vertex(2)->point();

      polygon.clear();
      polygon.push_back(Point_3(p0.x(), p0.y(), FT(0)));
      polygon.push_back(Point_3(p1.x(), p1.y(), FT(0)));
      polygon.push_back(Point_3(p2.x(), p2.y(), FT(0)));
      polygons.push_back(polygon);

      Color color(125, 125, 125);
      if (fit->info().index == KSR::uninitialized())
        color = Color(125, 125, 125);
      else if (fit->info().index == KSR::no_element())
        color = Color(125, 0, 0);
      else color = m_saver.get_idx_color(fit->info().index);
      colors.push_back(color);
    }
    m_saver.export_polygon_soup_3(polygons, colors, file_name);
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_POLYGON_SPLITTER_H
