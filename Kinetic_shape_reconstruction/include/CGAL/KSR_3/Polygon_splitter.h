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

// TODO: DOES NOT WORK WITH INEXACT KERNEL!
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
    std::vector<std::size_t> pedge_indices;
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
  using TAG = CGAL::Exact_intersections_tag;
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

public:
  Polygon_splitter(Data_structure& data) :
  m_data(data),
  m_merge_type(Planar_shape_type::CONVEX_HULL),
  m_verbose(m_data.is_verbose())
  { }

  void split_support_plane(const std::size_t sp_idx) {

    // if (sp_idx != 6) return;

    // Preprocessing.
    std::cout.precision(20);
    if (m_data.pfaces(sp_idx).size() > 1) merge_coplanar_pfaces(sp_idx);
    CGAL_assertion_msg(m_data.pfaces(sp_idx).size() == 1,
    "ERROR: WE CANNOT HAVE MULTIPLE COPLANAR PFACES!");
    const auto pface = *m_data.pfaces(sp_idx).begin();
    CGAL_assertion(pface.first == sp_idx);
    const auto original_input = m_data.input(pface);

    // Create cdt.
    initialize_cdt(pface);
    // dump_cdt(m_data, pface.first, m_cdt, "0-initial-");
    tag_cdt_exterior_faces();
    // dump_cdt(m_data, pface.first, m_cdt, "1-exterior-");
    tag_cdt_interior_faces();
    // dump_cdt(m_data, pface.first, m_cdt, "2-interior-");

    // Split polygons using cdt.
    m_data.clear_polygon_faces(sp_idx);
    initialize_new_pfaces(pface.first, original_input);

    // Set intersection adjacencies.
    reconnect_pvertices_to_ivertices();
    reconnect_pedges_to_iedges();
    set_new_adjacencies(pface.first);
  }

  void clear() {
    m_cdt.clear();
    m_input.clear();
    m_map_intersections.clear();
  }

private:
  Data_structure& m_data;
  TRI m_cdt;
  std::set<PVertex> m_input;
  std::map<CID, IEdge> m_map_intersections;
  const Planar_shape_type m_merge_type;
  const bool m_verbose;

  /*******************************
  **        MERGE PFACES        **
  ********************************/

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
    CGAL_assertion(is_pface_inside_bbox(support_plane_idx, merged));
  }

  // Check if the newly created pface goes beyond the bbox.
  const bool is_pface_inside_bbox(
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

    m_data.clear_pfaces(support_plane_idx);
    m_data.add_input_polygon(support_plane_idx, input_indices, merged);
  }

  /*******************************
  **         CREATE CDT         **
  ********************************/

  template<
  typename ForwardIt,
  typename BinaryPredicate>
  const ForwardIt unique_elements(
  ForwardIt first, ForwardIt last, const BinaryPredicate& predicate) const {

    if (first == last) return last;
    ForwardIt result = first;
    while (++first != last) {
      if (!predicate(*result, *first)) {
        if (++result != first) {
          *result = std::move(*first);
        }
      } else {
        auto& a = (*result).second;
        auto& b = (*first).second;
        if (
          a.first != Data_structure::null_pvertex() &&
          b.first == Data_structure::null_pvertex()) {
          b.first = a.first;
        }
        if (
          a.first == Data_structure::null_pvertex() &&
          b.first != Data_structure::null_pvertex()) {
          a.first = b.first;
        }
        if (
          a.second != Data_structure::null_ivertex() &&
          b.second == Data_structure::null_ivertex()) {
          b.second = a.second;
        }
        if (
          a.second == Data_structure::null_ivertex() &&
          b.second != Data_structure::null_ivertex()) {
          a.second = b.second;
        }
      }
    }
    return ++result;
  }

  void initialize_cdt(const PFace& pface) {

    std::cout.precision(20);
    const std::size_t sp_idx = pface.first;
    const auto& sp = m_data.support_plane(sp_idx);
    const auto pvertices = m_data.pvertices_of_pface(pface);
    const auto& iedges = sp.unique_iedges();

    // Create unique pvertices and ivertices.
    using Pair = std::pair<Point_2, std::pair<PVertex, IVertex> >;
    std::vector<Pair> points;
    points.reserve(pvertices.size() + iedges.size() * 2);

    for (const auto pvertex : pvertices) {
      const auto point = m_data.point_2(pvertex);
      points.push_back(std::make_pair(point,
        std::make_pair(pvertex, m_data.null_ivertex())));
    }

    for (const auto& iedge : iedges) {
      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);
      CGAL_assertion(isource != itarget);

      const auto source = m_data.to_2d(sp_idx, isource);
      const auto target = m_data.to_2d(sp_idx, itarget);
      CGAL_assertion(source != target);

      points.push_back(std::make_pair(source,
        std::make_pair(m_data.null_pvertex(), isource)));
      points.push_back(std::make_pair(target,
        std::make_pair(m_data.null_pvertex(), itarget)));
    }

    CGAL_assertion(points.size() == (pvertices.size() + iedges.size() * 2));
    // std::cout << "- num unique 1: " << points.size() << std::endl;

    const FT ptol = KSR::point_tolerance<FT>();
    const auto sort_cmp = [&](const Pair& a, const Pair& b) {
      const auto are_equal = ( KSR::distance(a.first, b.first) < ptol );
      if (!are_equal) return a.first < b.first;
      return false;
    };
    const auto unique_cmp = [&](const Pair& a, const Pair& b) {
      return ( KSR::distance(a.first, b.first) < ptol );
    };
    std::sort(points.begin(), points.end(), sort_cmp);
    points.erase(unique_elements(
      points.begin(), points.end(), unique_cmp), points.end());
    // std::cout << "- num unique 2: " << points.size() << std::endl;

    // for (const auto& pair : points) {
    //   std::cout <<
    //   m_data.str(pair.second.first) << " : " <<
    //   m_data.str(pair.second.second) << std::endl;
    //   std::cout << m_data.to_3d(sp_idx, pair.first) << std::endl;
    // }

    // Insert pvertices and ivertices.
    std::map<PVertex, Vertex_handle> vhs_pv;
    std::map<IVertex, Vertex_handle> vhs_iv;
    for (const auto& pair : points) {
      const auto& point = pair.first;
      const auto& data = pair.second;
      CGAL_assertion(
        data.first  != m_data.null_pvertex() ||
        data.second != m_data.null_ivertex());
      const auto vh = m_cdt.insert(point);

      if (data.first != m_data.null_pvertex()) {
        vh->info().pvertex = data.first;
        vhs_pv[vh->info().pvertex] = vh;
      }

      if (data.second != m_data.null_ivertex()) {
        vh->info().ivertex = data.second;
        vhs_iv[vh->info().ivertex] = vh;
      }
    }
    CGAL_assertion(vhs_pv.size() >= 3);
    CGAL_assertion(vhs_iv.size() >= 1);
    // std::cout << "- num cdt verts 1: " << m_cdt.number_of_vertices() << std::endl;
    // std::cout << "- num cdt faces 1: " << m_cdt.number_of_faces()    << std::endl;

    // Insert pedge constraints.
    std::vector<PVertex> polygon;
    polygon.reserve(pvertices.size());
    std::copy(pvertices.begin(), pvertices.end(), std::back_inserter(polygon));
    CGAL_assertion(polygon.size() == pvertices.size());
    const std::size_t n = polygon.size();
    std::vector< std::vector<IVertex> > pedge_map(n);

    for (const auto& pair : points) {
      const auto& point = pair.first;
      const auto& data = pair.second;
      if (data.first != m_data.null_pvertex()) continue;
      CGAL_assertion(data.first == m_data.null_pvertex());
      CGAL_assertion(data.second != m_data.null_ivertex());
      CGAL_assertion(!is_pvertex(vhs_pv, polygon, point));

      const std::size_t idx = find_pedge(vhs_pv, polygon, point);
      if (idx != KSR::no_element()) {
        CGAL_assertion(idx < pedge_map.size());
        pedge_map[idx].push_back(data.second);
      }
    }

    CGAL_assertion(pedge_map.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
      if (pedge_map[i].size() > 0) continue;
      const std::size_t ip = (i + 1) % n;
      const auto& psource = polygon[i];
      const auto& ptarget = polygon[ip];

      CGAL_assertion(psource != ptarget);
      CGAL_assertion(vhs_pv.find(psource) != vhs_pv.end());
      CGAL_assertion(vhs_pv.find(ptarget) != vhs_pv.end());

      const auto& vh_source = vhs_pv.at(psource);
      const auto& vh_target = vhs_pv.at(ptarget);
      CGAL_assertion(vh_source != vh_target);
      CGAL_assertion(KSR::distance(
        vh_source->point(), vh_target->point()) >= ptol);
      m_cdt.insert_constraint(vh_source, vh_target);
      m_input.insert(psource);
    }

    // Set pedge indices.
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t im = (i + n - 1) % n;
      const auto& pvertex = polygon[i];
      CGAL_assertion(vhs_pv.find(pvertex) != vhs_pv.end());
      const auto& vh = vhs_pv.at(pvertex);
      CGAL_assertion(vh->info().pedge_indices.size() == 0);
      vh->info().pedge_indices.push_back(im);
      vh->info().pedge_indices.push_back(i);
    }

    for (std::size_t i = 0; i < pedge_map.size(); ++i) {
      const auto& ivertices = pedge_map[i];
      if (ivertices.size() == 0) continue;
      for (const auto& ivertex : ivertices) {
        CGAL_assertion(vhs_iv.find(ivertex) != vhs_iv.end());
        const auto& vh = vhs_iv.at(ivertex);
        if (vh->info().pvertex != m_data.null_pvertex()) {
          CGAL_assertion(vh->info().pedge_indices.size() == 2);
          continue;
        }
        CGAL_assertion(vh->info().pvertex == m_data.null_pvertex());
        CGAL_assertion(vh->info().pedge_indices.size() == 0);
        vh->info().pedge_indices.push_back(i);
      }
    }

    // std::cout << "- num cdt verts 2: " << m_cdt.number_of_vertices() << std::endl;
    // std::cout << "- num cdt faces 2: " << m_cdt.number_of_faces()    << std::endl;

    // Insert iedge constraints.
    for (const auto& iedge : iedges) {
      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);

      CGAL_assertion(isource != itarget);
      CGAL_assertion(vhs_iv.find(isource) != vhs_iv.end());
      CGAL_assertion(vhs_iv.find(itarget) != vhs_iv.end());

      const auto& vh_source = vhs_iv.at(isource);
      const auto& vh_target = vhs_iv.at(itarget);
      CGAL_assertion(vh_source != vh_target);
      CGAL_assertion(KSR::distance(
        vh_source->point(), vh_target->point()) >= ptol);

      const auto cid = m_cdt.insert_constraint(vh_source, vh_target);
      CGAL_assertion(m_map_intersections.find(cid) == m_map_intersections.end());
      m_map_intersections.insert(std::make_pair(cid, iedge));
    }

    // std::cout << "- num cdt verts 3: " << m_cdt.number_of_vertices() << std::endl;
    // std::cout << "- num cdt faces 3: " << m_cdt.number_of_faces()    << std::endl;
  }

  const bool is_pvertex(
    const std::map<PVertex, Vertex_handle>& vhs_pv,
    const std::vector<PVertex>& polygon,
    const Point_2& query) const {

    const FT ptol = KSR::point_tolerance<FT>();
    for (const auto& pvertex : polygon) {
      CGAL_assertion(vhs_pv.find(pvertex) != vhs_pv.end());
      const auto& vh = vhs_pv.at(pvertex);
      const auto& point = vh->point();
      const FT distance = KSR::distance(point, query);
      if (distance < ptol) return true;
    }
    return false;
  }

  const std::size_t find_pedge(
    const std::map<PVertex, Vertex_handle>& vhs_pv,
    const std::vector<PVertex>& polygon,
    const Point_2& query) const {

    const FT tol = KSR::tolerance<FT>();
    const FT ptol = KSR::point_tolerance<FT>();

    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& psource = polygon[i];
      const auto& ptarget = polygon[ip];

      CGAL_assertion(psource != ptarget);
      CGAL_assertion(vhs_pv.find(psource) != vhs_pv.end());
      CGAL_assertion(vhs_pv.find(ptarget) != vhs_pv.end());

      const auto& vh_source = vhs_pv.at(psource);
      const auto& vh_target = vhs_pv.at(ptarget);
      CGAL_assertion(vh_source != vh_target);

      const auto& source = vh_source->point();
      const auto& target = vh_target->point();
      CGAL_assertion(KSR::distance(source, target) >= ptol);

      const FT half = FT(1) / FT(2);
      const Vector_2 s1(query, source);
      const Vector_2 s2(query, target);

      const FT A = half * CGAL::determinant(s1, s2);
      const FT D = CGAL::scalar_product(s1, s2);
      if (CGAL::abs(A) < tol && D < FT(0)) return i;
    }
    return KSR::no_element();
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
        const bool is_boundary_edge = is_boundary(edge);
        if (!is_boundary_edge) {
          todo.push(next);
        }
      }
    }
    CGAL_assertion(todo.size() == 0);
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

  const bool is_boundary(const Edge& edge) const {
    const auto& fh = edge.first;
    const std::size_t idx = edge.second;

    const auto& vh1 = fh->vertex( (idx + 1) % 3 );
    const auto& vh2 = fh->vertex( (idx + 2) % 3 );

    const auto& pes1 = vh1->info().pedge_indices;
    const auto& pes2 = vh2->info().pedge_indices;
    CGAL_assertion(pes1.size() <= 2);
    CGAL_assertion(pes2.size() <= 2);

    if (pes1.size() == 0) return false;
    if (pes2.size() == 0) return false;
    CGAL_assertion(pes1.size() > 0);
    CGAL_assertion(pes2.size() > 0);

    for (const std::size_t pe1 : pes1) {
      for (const std::size_t pe2 : pes2) {
        if (pe1 == pe2) return true;
      }
    }
    return false;
  }

  void initialize_new_pfaces(
    const std::size_t sp_idx, const std::vector<std::size_t>& original_input) {

    std::size_t num_pfaces = 0;
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
        if (source->info().pvertex == m_data.null_pvertex()) {
          const auto& p = source->point();
          const Point_2 spoint(
            static_cast<FT>(CGAL::to_double(p.x())),
            static_cast<FT>(CGAL::to_double(p.y())));
          source->info().pvertex = m_data.add_pvertex(sp_idx, spoint);
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
      const auto pface = m_data.add_pface(new_pvertices);
      CGAL_assertion(pface != PFace());
      m_data.input(pface) = original_input;
      ++num_pfaces;
    }

    if (m_verbose) {
      std::cout << "- number of newly inserted pfaces: " << num_pfaces << std::endl;
    }
  }

  void reconnect_pvertices_to_ivertices() {

    // Reconnect only those, which have already been connected.
    for (auto vit = m_cdt.finite_vertices_begin(); vit != m_cdt.finite_vertices_end(); ++vit) {
      if (vit->info().pvertex != m_data.null_pvertex() &&
          vit->info().ivertex != m_data.null_ivertex()) {
        m_data.connect(vit->info().pvertex, vit->info().ivertex);
      }
    }
  }

  void reconnect_pedges_to_iedges() {

    // Reconnect only those, which have already been connected.
    for (const auto& item : m_map_intersections) {
      const auto& cid   = item.first;
      const auto& iedge = item.second;

      if (iedge == m_data.null_iedge()) continue;
      CGAL_assertion(iedge != m_data.null_iedge());

      auto vit = m_cdt.vertices_in_constraint_begin(cid);
      while (true) {
        auto next = vit; ++next;
        if (next == m_cdt.vertices_in_constraint_end(cid)) { break; }
        const auto a = *vit;
        const auto b = *next;
        vit = next;

        if (
          a->info().pvertex == m_data.null_pvertex() ||
          b->info().pvertex == m_data.null_pvertex()) {
          continue;
        }
        CGAL_assertion(a->info().pvertex != m_data.null_pvertex());
        CGAL_assertion(b->info().pvertex != m_data.null_pvertex());
        m_data.connect(a->info().pvertex, b->info().pvertex, iedge);
      }
    }
  }

  void set_new_adjacencies(const std::size_t sp_idx) {

    // std::cout << std::endl << "support plane idx: " << sp_idx << std::endl;
    const auto all_pvertices = m_data.pvertices(sp_idx);
    for (const auto pvertex : all_pvertices) {
      // std::cout << "pvertex: " << m_data.point_3(pvertex) << std::endl;

      bool is_frozen = false;
      auto iedge = m_data.null_iedge();
      std::pair<PVertex, PVertex> neighbors(
        m_data.null_pvertex(), m_data.null_pvertex());

      // Search for a frozen pvertex.
      const auto pedges = m_data.pedges_around_pvertex(pvertex);
      for (const auto pedge : pedges) {
        // std::cout << "pedge: 2 " << m_data.segment_3(pedge) << " : "
        // << m_data.has_iedge(pedge) << std::endl;

        if (m_data.has_iedge(pedge)) {
          if (iedge == m_data.null_iedge()) {
            // std::cout << "empty iedge" << std::endl;
            iedge = m_data.iedge(pedge);
          } else {
            // std::cout << "frozen pvertex" << std::endl;
            is_frozen = true;
            break;
          }
        } else {
          const auto opposite = m_data.opposite(pedge, pvertex);
          if (neighbors.first == m_data.null_pvertex()) {
            neighbors.first = opposite;
            // std::cout << "assigned first neighbor: " << m_data.point_3(opposite) << std::endl;
          } else {
            CGAL_assertion(neighbors.first  != m_data.null_pvertex());
            CGAL_assertion(neighbors.second == m_data.null_pvertex());
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
      if (iedge == m_data.null_iedge()) {
        continue;
      }
      m_data.connect(pvertex, iedge);
      // CGAL_assertion(
      //   neighbors.first  != m_data.null_pvertex() &&
      //   neighbors.second != m_data.null_pvertex());

      // Set future direction.
      bool is_first_okay = false;
      if (neighbors.first != m_data.null_pvertex()) {
        is_first_okay = update_neighbor(pvertex,   neighbors.first);
      }

      bool is_second_okay = false;
      if (neighbors.second != m_data.null_pvertex()) {
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

      const auto intersection_line = m_data.segment_2(sp_idx, iedge).supporting_line();
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
};

} // namespace KSR_3
} // namespace CGAL

#endif // CGAL_KSR_3_POLYGON_SPLITTER_H
