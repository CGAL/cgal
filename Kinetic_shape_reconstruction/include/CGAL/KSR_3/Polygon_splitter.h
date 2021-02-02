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
#include <CGAL/convex_hull_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/KSR/enum.h>
#include <CGAL/KSR/utils.h>
#include <CGAL/KSR/debug.h>

namespace CGAL {
namespace KSR_3 {

template<typename Data_structure, typename GeomTraits>
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
  using Segment_2  = typename Kernel::Segment_2;

  using PVertex = typename Data_structure::PVertex;
  using PFace   = typename Data_structure::PFace;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  struct Vertex_info {
    PVertex pvertex;
    IVertex ivertex;
    Vertex_info() :
    pvertex(Data_structure::null_pvertex()),
    ivertex(Data_structure::null_ivertex())
    { }
  };

  struct Face_info {
    bool tagged;
    std::size_t index;
    std::size_t input;
    Face_info() :
    tagged(false),
    index(KSR::uninitialized()),
    input(KSR::uninitialized())
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

  using Planar_shape_type = KSR::Planar_shape_type;

  Data_structure& m_data;
  TRI m_cdt;
  std::set<PVertex> m_input;
  std::map<CID, IEdge> m_map_intersections;
  const Planar_shape_type m_merge_type;

public:
  Polygon_splitter(Data_structure& data) :
  m_data(data),
  m_merge_type(Planar_shape_type::CONVEX_HULL)
  { }

  void split_support_plane(const std::size_t support_plane_idx) {

    // Preprocessing.
    const auto all_pfaces = m_data.pfaces(support_plane_idx);
    if (all_pfaces.size() > 1) {
      merge_coplanar_pfaces(support_plane_idx);
    }

    // Create cdt.
    std::cout.precision(20);
    std::vector< std::vector<std::size_t> > original_input;
    std::vector< std::vector<Point_2> > original_faces;
    initialize_cdt(support_plane_idx, original_input, original_faces);
    CGAL_assertion(original_faces.size() >= 1);
    CGAL_assertion(original_input.size() == original_faces.size());
    tag_cdt_exterior_faces();
    tag_cdt_interior_faces();

    // Split polygons using cdt.
    m_data.clear_polygon_faces(support_plane_idx);
    initialize_new_pfaces(
      support_plane_idx, original_input, original_faces);

    // Set intersection adjacencies.
    reconnect_pvertices_to_ivertices();
    reconnect_pedges_to_iedges();
    set_new_adjacencies(support_plane_idx);
    if (original_faces.size() > 1) {
      CGAL_assertion_msg(false,
      "ERROR: WE CANNOT HAVE MULTIPLE COPLANAR PFACES!");
    }
  }

private:

  void merge_coplanar_pfaces(
    const std::size_t support_plane_idx) {

    const bool is_debug = false;
    CGAL_assertion(support_plane_idx >= 6);
    if (is_debug) {
      std::cout << std::endl << "support plane idx: " << support_plane_idx << std::endl;
    }

    std::vector<Point_2> points;
    collect_pface_points(support_plane_idx, points);
    std::vector<Point_2> merged;
    create_merged_pface(support_plane_idx, points, merged);

    if (is_debug) {
      std::cout << "merged pface: " << std::endl;
      for (std::size_t i = 0; i < merged.size(); ++i) {
        const std::size_t ip = (i + 1) % merged.size();
        const auto& p = merged[i];
        const auto& q = merged[ip];
        std::cout << "2 " <<
        m_data.to_3d(support_plane_idx, p) << " " <<
        m_data.to_3d(support_plane_idx, q) << std::endl;
      }
    }
    add_merged_pface(support_plane_idx, merged);
  }

  void collect_pface_points(
    const std::size_t support_plane_idx,
    std::vector<Point_2>& points) const {

    points.clear();
    const auto all_pfaces = m_data.pfaces(support_plane_idx);
    CGAL_assertion(all_pfaces.size() > 1);
    points.reserve(all_pfaces.size() * 3);
    for (const auto pface : all_pfaces) {
      const auto pvertices = m_data.pvertices_of_pface(pface);
      CGAL_assertion(pvertices.size() >= 3);

      for (const auto pvertex : pvertices) {
        const auto point = m_data.point_2(pvertex);
        points.push_back(point);
      }
    }
    CGAL_assertion(points.size() >= all_pfaces.size() * 3);
  }

  void create_merged_pface(
    const std::size_t support_plane_idx,
    const std::vector<Point_2>& points,
    std::vector<Point_2>& merged) const {

    merged.clear();
    switch (m_merge_type) {
      case Planar_shape_type::CONVEX_HULL: {
        CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(merged) );
        break;
      }
      case Planar_shape_type::RECTANGLE: {
        CGAL_assertion_msg(false, "TODO: MERGE PFACES INTO A RECTANGLE!");
        break;
      }
      default: {
        CGAL_assertion_msg(false, "ERROR: MERGE PFACES, WRONG TYPE!");
        break;
      }
    }
    CGAL_assertion(merged.size() >= 3);
    CGAL_assertion(check_merged_pface(support_plane_idx, merged));
  }

  // Check if the newly created pface goes beyond the bbox.
  const bool check_merged_pface(
    const std::size_t support_plane_idx,
    const std::vector<Point_2>& merged) const {

    std::vector<Point_2> bbox;
    create_bbox(support_plane_idx, bbox);
    CGAL_assertion(bbox.size() == 4);

    for (std::size_t i = 0; i < 4; ++i) {
      const std::size_t ip = (i + 1) % 4;
      const auto& pi = bbox[i];
      const auto& qi = bbox[ip];
      const Segment_2 edge(pi, qi);

      for (std::size_t j = 0; j < merged.size(); ++j) {
        const std::size_t jp = (j + 1) % merged.size();
        const auto& pj = merged[j];
        const auto& qj = merged[jp];
        const Segment_2 segment(pj, qj);
        Point_2 inter;
        const bool is_intersected = KSR::intersection(segment, edge, inter);
        if (is_intersected) return false;
      }
    }
    return true;
  }

  void add_merged_pface(
    const std::size_t support_plane_idx,
    std::vector<Point_2>& merged) {

    const auto all_pfaces = m_data.pfaces(support_plane_idx);
    std::vector<std::size_t> input_indices;
    input_indices.reserve(all_pfaces.size());

    for (const auto pface : all_pfaces) {
      const auto& pface_input = m_data.input(pface);
      CGAL_assertion(pface_input.size() == 1);
      input_indices.push_back(pface_input[0]);
    }
    CGAL_assertion(input_indices.size() == all_pfaces.size());

    m_data.support_plane(support_plane_idx).clear_pfaces();
    m_data.add_input_polygon(support_plane_idx, input_indices, merged);
  }

  void initialize_cdt(
    const std::size_t support_plane_idx,
    std::vector< std::vector<std::size_t> >& original_input,
    std::vector< std::vector<Point_2> >& original_faces) {

    // Insert pvertices.
    std::map<PVertex, Vertex_handle> vhs_map;
    const auto all_pvertices = m_data.pvertices(support_plane_idx);
    for (const auto pvertex : all_pvertices) {
      const auto vh = m_cdt.insert(m_data.point_2(pvertex));
      vh->info().pvertex = pvertex;
      m_input.insert(pvertex);
      vhs_map[pvertex] = vh;
    }

    // Insert pfaces and the corresponding constraints.
    original_faces.clear();
    original_input.clear();

    // std::vector<Vertex_handle> vhs;
    std::vector<Point_2> original_face;
    const auto all_pfaces = m_data.pfaces(support_plane_idx);
    for (const auto pface : all_pfaces) {
      const auto pvertices = m_data.pvertices_of_pface(pface);

      // vhs.clear();
      original_face.clear();
      for (const auto pvertex : pvertices) {
        CGAL_assertion(vhs_map.find(pvertex) != vhs_map.end());
        const auto vh = vhs_map.at(pvertex);
        original_face.push_back(vh->point());
        // vhs.push_back(vh);
      }

      original_faces.push_back(original_face);
      original_input.push_back(m_data.input(pface));

      // TODO: WHY WE CANNOT USE VHS DIRECTLY HERE? THAT SHOULD BE MORE PRECISE!
      // vhs.push_back(vhs.front());
      // const auto cid = m_cdt.insert_constraint(vhs.begin(), vhs.end());
      original_face.push_back(original_face.front());
      const auto cid = m_cdt.insert_constraint(original_face.begin(), original_face.end());
      m_map_intersections.insert(std::make_pair(cid, Data_structure::null_iedge()));
    }

    // Then, add intersection vertices + constraints.
    const auto& iedges = m_data.support_plane(support_plane_idx).unique_iedges();
    for (const auto& iedge : iedges) {
      const auto source = m_data.source(iedge);
      const auto target = m_data.target(iedge);

      const auto vsource = m_cdt.insert(m_data.to_2d(support_plane_idx, source));
      vsource->info().ivertex = source;
      const auto vtarget = m_cdt.insert(m_data.to_2d(support_plane_idx, target));
      vtarget->info().ivertex = target;

      const auto cid = m_cdt.insert_constraint(vsource, vtarget);
      m_map_intersections.insert(std::make_pair(cid, iedge));
    }

    // Finally, add original labels to the cdt.
    if (all_pfaces.size() > 1) {
      CGAL_assertion_msg(false,
      "ERROR: WE CANNOT HAVE MULTIPLE COPLANAR PFACES!");
    }
  }

  // All exterior faces are tagged by KSR::no_element().
  void tag_cdt_exterior_faces() {

    std::queue<Face_handle> todo;
    todo.push(m_cdt.incident_faces(m_cdt.infinite_vertex()));
    while (!todo.empty()) {
      const auto fh = todo.front();
      todo.pop();
      if (fh->info().index != KSR::uninitialized()) {
        continue;
      }
      fh->info().index = KSR::no_element();

      for (std::size_t i = 0; i < 3; ++i) {
        const auto next = fh->neighbor(i);
        const auto edge = std::make_pair(fh, i);
        const bool is_border_edge = is_border(edge);
        if (!is_border_edge) {
          todo.push(next);
        }
      }
    }
    CGAL_assertion(todo.size() == 0);
  }

  const bool is_border(const Edge& edge) const {

    if (!m_cdt.is_constrained(edge)) {
      return false;
    }

    const std::size_t im = (edge.second + 2) % 3;
    const std::size_t ip = (edge.second + 1) % 3;

    const auto vm = edge.first->vertex(im);
    const auto vp = edge.first->vertex(ip);

    const auto ctx_begin = m_cdt.contexts_begin(vp, vm);
    const auto ctx_end = m_cdt.contexts_end(vp, vm);

    for (auto cit = ctx_begin; cit != ctx_end; ++cit) {
      const auto iter = m_map_intersections.find(cit->id());
      if (iter == m_map_intersections.end()) {
        continue;
      }
      if (iter->second == Data_structure::null_iedge()) {
        return true;
      }
    }
    return false;
  }

  // All enterior faces are tagged by face_index.
  void tag_cdt_interior_faces() {

    std::size_t face_index = 0;
    std::queue<Face_handle> todo;
    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      CGAL_assertion(todo.size() == 0);
      if (fit->info().index != KSR::uninitialized()) {
        continue;
      }

      todo.push(fit);
      std::size_t num_faces = 0;
      while (!todo.empty()) {
        const auto fh = todo.front();
        todo.pop();
        if (fh->info().index != KSR::uninitialized()) {
          continue;
        }
        fh->info().index = face_index;
        ++num_faces;

        for (std::size_t i = 0; i < 3; ++i) {
          const auto next = fh->neighbor(i);
          const auto edge = std::make_pair(fh, i);
          const bool is_constrained_edge = m_cdt.is_constrained(edge);
          if (!is_constrained_edge) {
            todo.push(next);
          }
        }
      }
      ++face_index;
      CGAL_assertion(todo.size() == 0);
    }
  }

  void initialize_new_pfaces(
    const std::size_t support_plane_idx,
    const std::vector< std::vector<std::size_t> >& original_input,
    const std::vector< std::vector<Point_2> >& original_faces) {

    std::set<std::size_t> done;
    for (auto fit = m_cdt.finite_faces_begin(); fit != m_cdt.finite_faces_end(); ++fit) {
      CGAL_assertion(fit->info().index != KSR::uninitialized());
      if (fit->info().index == KSR::no_element()) { // skip exterior faces
        continue;
      }

      // Search for a constrained edge.
      Edge edge;
      for (std::size_t i = 0; i < 3; ++i) {
        edge = std::make_pair(fit, i);
        if (m_cdt.is_constrained(edge)) {
          break;
        }
      }

      // Skip pure interior faces.
      if (!m_cdt.is_constrained(edge)) {
        continue;
      }

      // If face index is already a part of the set, skip.
      const auto fh = edge.first;
      if (!done.insert(fh->info().index).second) {
        continue;
      }

      // Start from the constrained edge and traverse all constrained edges / boundary
      // of the triangulation part that is tagged with the same face index.
      // While traversing, add all missing pvertices.
      auto curr = edge;
      std::vector<PVertex> new_pvertices;
      do {
        const auto curr_face = curr.first;
        const int idx = curr.second;

        const auto source = curr_face->vertex(m_cdt.ccw(idx));
        const auto target = curr_face->vertex(m_cdt.cw (idx));
        if (source->info().pvertex == Data_structure::null_pvertex()) {
          source->info().pvertex =
            m_data.add_pvertex(support_plane_idx, source->point());
        }
        new_pvertices.push_back(source->info().pvertex);

        // Search for the next constrained edge.
        auto next = std::make_pair(curr_face, m_cdt.ccw(idx));
        while (!m_cdt.is_constrained(next)) {

          const auto next_face = next.first->neighbor(next.second);
          // Should be the same original polygon.
          CGAL_assertion(next_face->info().index == edge.first->info().index);

          const int next_idx = m_cdt.ccw(next_face->index(next.first));
          next = std::make_pair(next_face, next_idx);
        }
        // Check wether next source == previous target.
        CGAL_assertion(next.first->vertex(m_cdt.ccw(next.second)) == target);
        curr = next;

      } while (curr != edge);
      CGAL_assertion(curr == edge);

      // Add a new pface.
      CGAL_assertion(original_faces.size() > 0);
      CGAL_assertion(original_input.size() == original_faces.size());
      const auto pface = m_data.add_pface(new_pvertices);
      CGAL_assertion(pface != PFace());

      if (original_faces.size() != 1) {
        CGAL_assertion_msg(original_faces.size() <= 1,
        "ERROR: WE CANNOT HAVE MULTIPLE COPLANAR PFACES!");
      } else {
        m_data.input(pface) = original_input[0];
      }
    }
  }

  void reconnect_pvertices_to_ivertices() {

    // Reconnect only those, which have already been connected.
    for (auto vit = m_cdt.finite_vertices_begin(); vit != m_cdt.finite_vertices_end(); ++vit) {
      if (vit->info().pvertex != Data_structure::null_pvertex() &&
          vit->info().ivertex != Data_structure::null_ivertex()) {
        m_data.connect(vit->info().pvertex, vit->info().ivertex);
      }
    }
  }

  void reconnect_pedges_to_iedges() {

    // Reconnect only those, which have already been connected.
    for (const auto& item : m_map_intersections) {
      const auto& cid   = item.first;
      const auto& iedge = item.second;

      if (iedge == Data_structure::null_iedge()) {
        continue;
      }
      CGAL_assertion(iedge != Data_structure::null_iedge());

      auto vit = m_cdt.vertices_in_constraint_begin(cid);
      while (true) {
        auto next = vit; ++next;
        if (next == m_cdt.vertices_in_constraint_end(cid)) { break; }
        const auto a = *vit;
        const auto b = *next;
        vit = next;

        if (
          a->info().pvertex == Data_structure::null_pvertex() ||
          b->info().pvertex == Data_structure::null_pvertex()) {
          continue;
        }
        CGAL_assertion(a->info().pvertex != Data_structure::null_pvertex());
        CGAL_assertion(b->info().pvertex != Data_structure::null_pvertex());
        m_data.connect(a->info().pvertex, b->info().pvertex, iedge);
      }
    }
  }

  void set_new_adjacencies(
    const std::size_t support_plane_idx) {

    // std::cout << std::endl << "support plane idx: " << support_plane_idx << std::endl;
    const auto all_pvertices = m_data.pvertices(support_plane_idx);
    for (const auto pvertex : all_pvertices) {
      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;

      bool is_frozen = false;
      auto iedge = Data_structure::null_iedge();
      std::pair<PVertex, PVertex> neighbors(
        Data_structure::null_pvertex(), Data_structure::null_pvertex());

      // Search for a frozen pvertex.
      const auto pedges = m_data.pedges_around_pvertex(pvertex);
      for (const auto pedge : pedges) {
        // std::cout << "pedge: 2 " << m_data.segment_3(pedge) << " : "
        // << m_data.has_iedge(pedge) << std::endl;

        if (m_data.has_iedge(pedge)) {
          if (iedge == Data_structure::null_iedge()) {
            // std::cout << "empty iedge" << std::endl;
            iedge = m_data.iedge(pedge);
          } else {
            // std::cout << "frozen pvertex" << std::endl;
            is_frozen = true;
            break;
          }
        } else {
          const auto opposite = m_data.opposite(pedge, pvertex);
          if (neighbors.first == Data_structure::null_pvertex()) {
            neighbors.first = opposite;
            // std::cout << "assigned first neighbor: " << m_data.point_3(opposite) << std::endl;
          } else {
            CGAL_assertion(neighbors.first  != Data_structure::null_pvertex());
            CGAL_assertion(neighbors.second == Data_structure::null_pvertex());
            neighbors.second = opposite;
            // std::cout << "assigned second neighbor: " << m_data.point_3(opposite) << std::endl;
          }
        }
      }

      // Several incident intersections = frozen pvertex.
      if (is_frozen) {
        m_data.direction(pvertex) = CGAL::NULL_VECTOR;
        continue;
      }

      // No incident intersections = keep initial direction.
      if (iedge == Data_structure::null_iedge()) {
        continue;
      }
      m_data.connect(pvertex, iedge);
      // CGAL_assertion(
      //   neighbors.first  != Data_structure::null_pvertex() &&
      //   neighbors.second != Data_structure::null_pvertex());

      // Set future direction.
      bool is_first_okay = false;
      if (neighbors.first != Data_structure::null_pvertex()) {
        is_first_okay = update_neighbor(pvertex,   neighbors.first);
      }

      bool is_second_okay = false;
      if (neighbors.second != Data_structure::null_pvertex()) {
        is_second_okay = update_neighbor(pvertex, neighbors.second);
      }

      Line_2 future_line;
      if (is_first_okay && is_second_okay) {
        future_line = Line_2(
          m_data.point_2(neighbors.first , FT(1)),
          m_data.point_2(neighbors.second, FT(1)));
      } else {
        CGAL_assertion(is_first_okay && !is_second_okay);
        future_line = Line_2(
          m_data.point_2(pvertex        , FT(1)),
          m_data.point_2(neighbors.first, FT(1)));
      }
      CGAL_assertion(future_line != Line_2());

      const auto intersection_line = m_data.segment_2(support_plane_idx, iedge).supporting_line();
      CGAL_assertion_msg(!CGAL::parallel(intersection_line, future_line),
      "TODO: POLYGON SPLITTER, HANDLE CASE WITH PARALLEL LINES!");
      const Point_2 future_point = KSR::intersection<Point_2>(intersection_line, future_line);
      const auto pinit = m_data.point_2(pvertex, FT(0));
      const Vector_2 future_direction(pinit, future_point);
      m_data.direction(pvertex) = future_direction;
      // std::cout << "future point: " << m_data.to_3d(pvertex.first, future_point) << std::endl;
    }
  }

  const bool update_neighbor(
    const PVertex& pvertex, PVertex& neighbor) const {

    bool is_okay = (m_input.find(neighbor) != m_input.end());
    auto last = pvertex;
    auto curr = neighbor;
    while (!is_okay) {
      PVertex next, ignored;
      std::tie(next, ignored) = m_data.border_prev_and_next(curr);
      if (next == last) {
        std::swap(next, ignored);
      }
      CGAL_assertion(ignored == last);

      last = curr; curr = next;
      if (m_input.find(curr) != m_input.end()) {
        neighbor = curr;
        is_okay = true;
      }
    }
    return is_okay;
  }

  void create_bbox(
    const std::size_t support_plane_idx,
    std::vector<Point_2>& bbox) const {

    CGAL_assertion(support_plane_idx >= 6);
    const auto& iedges = m_data.support_plane(support_plane_idx).unique_iedges();

    std::vector<Point_2> points;
    points.reserve(iedges.size() * 2);

    for (const auto& iedge : iedges) {
      const auto source = m_data.source(iedge);
      const auto target = m_data.target(iedge);
      points.push_back(m_data.to_2d(support_plane_idx, source));
      points.push_back(m_data.to_2d(support_plane_idx, target));
    }
    CGAL_assertion(points.size() == iedges.size() * 2);

    const auto box = CGAL::bbox_2(points.begin(), points.end());
    const Point_2 p1(box.xmin(), box.ymin());
    const Point_2 p2(box.xmax(), box.ymin());
    const Point_2 p3(box.xmax(), box.ymax());
    const Point_2 p4(box.xmin(), box.ymax());

    bbox.clear();
    bbox.reserve(4);
    bbox.push_back(p1);
    bbox.push_back(p2);
    bbox.push_back(p3);
    bbox.push_back(p4);
  }

  void dump(
    const bool dump_data,
    const std::size_t type, // 0 - index, 1 - input
    const std::size_t support_plane_idx,
    std::string file_name = "") {
    if (!dump_data) return;

    Mesh_3 mesh;
    Uchar_map red   = mesh.template add_property_map<Face_index, unsigned char>("red"  , 125).first;
    Uchar_map green = mesh.template add_property_map<Face_index, unsigned char>("green", 125).first;
    Uchar_map blue  = mesh.template add_property_map<Face_index, unsigned char>("blue" , 125).first;

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
      if (type == 0) {
        CGAL::Random rand(fit->info().index);
        if (fit->info().index != KSR::no_element()) {
          red[face]   = (unsigned char)(rand.get_int(32, 192));
          green[face] = (unsigned char)(rand.get_int(32, 192));
          blue[face]  = (unsigned char)(rand.get_int(32, 192));
        }
      } else if (type == 1) {
        CGAL::Random rand(fit->info().input);
        if (fit->info().input != KSR::uninitialized()) {
          red[face]   = (unsigned char)(rand.get_int(32, 192));
          green[face] = (unsigned char)(rand.get_int(32, 192));
          blue[face]  = (unsigned char)(rand.get_int(32, 192));
        }
      } else {
        CGAL_assertion_msg(false, "ERROR: WRONG LABEL TYPE!");
      }
    }

    if (file_name == "")
      file_name = "support_cdt_" + std::to_string(support_plane_idx) + ".ply";
    std::ofstream out(file_name);
    out.precision(20);
    CGAL::write_ply(out, mesh);
    out.close();
  }

  void dump_original_faces(
    const std::size_t support_plane_idx,
    const std::vector< std::vector<Point_2> >& original_faces,
    const std::string file_name) {

    std::vector<Point_3> polygon;
    std::vector< std::vector<Point_3> > polygons;
    polygons.reserve(original_faces.size());

    for (const auto& original_face : original_faces) {
      polygon.clear();
      for (const auto& p : original_face) {
        polygon.push_back(m_data.to_3d(support_plane_idx, p));
      }
      polygons.push_back(polygon);
    }
    CGAL_assertion(polygons.size() == original_faces.size());
    KSR_3::Saver<Kernel> saver;
    saver.export_polygon_soup_3(polygons, file_name);
  }

  void dump_current_pface(
    const PFace& pface,
    const std::string file_name) {

    const auto pvertices = m_data.pvertices_of_pface(pface);
    std::vector< std::vector<Point_3> > polygons(1);
    for (const auto pvertex : pvertices)
      polygons[0].push_back(m_data.point_3(pvertex, FT(0)));
    KSR_3::Saver<Kernel> saver;
    saver.export_polygon_soup_3(polygons, file_name);
  }
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_POLYGON_SPLITTER_H
