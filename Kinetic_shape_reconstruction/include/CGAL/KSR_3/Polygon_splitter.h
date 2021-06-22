// Copyright (c) 2019 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
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
  using FT          = typename Kernel::FT;
  using Point_2     = typename Kernel::Point_2;
  using Point_3     = typename Kernel::Point_3;
  using Line_2      = typename Kernel::Line_2;
  using Vector_2    = typename Kernel::Vector_2;
  using Triangle_2  = typename Kernel::Triangle_2;
  using Segment_2   = typename Kernel::Segment_2;
  using Direction_2 = typename Kernel::Direction_2;

  using PVertex = typename Data_structure::PVertex;
  using PFace   = typename Data_structure::PFace;
  using PEdge   = typename Data_structure::PEdge;

  using IVertex = typename Data_structure::IVertex;
  using IEdge   = typename Data_structure::IEdge;

  struct Vertex_info {
    PVertex pvertex;
    IVertex ivertex;
    Vertex_info() :
    pvertex(Data_structure::null_pvertex()),
    ivertex(Data_structure::null_ivertex()),
    sp_idx(KSR::no_element()),
    is_boundary_ivertex(false)
    { }
    std::size_t sp_idx;
    std::vector<std::size_t> pedge_indices;
    bool is_boundary_ivertex;
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

    // if (sp_idx != 17) return;

    // Preprocessing.
    std::cout.precision(20);
    if (m_data.pfaces(sp_idx).size() > 1) {
      CGAL_assertion_msg(false, "ERROR: THIS CALL SHOULD NEVER HAPPEN!");
      merge_coplanar_pfaces(sp_idx);
    }
    CGAL_assertion_msg(m_data.pfaces(sp_idx).size() == 1,
    "ERROR: WE CANNOT HAVE MULTIPLE COPLANAR PFACES!");
    const auto pface = *m_data.pfaces(sp_idx).begin();
    CGAL_assertion(pface.first == sp_idx);
    const auto original_input = m_data.input(pface);
    CGAL_assertion(m_data.pvertices_of_pface(pface).size() >= 3);

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
    m_boundary_ivertices.clear();
  }

private:
  Data_structure& m_data;
  TRI m_cdt;
  std::set<PVertex> m_input;
  std::map<CID, IEdge> m_map_intersections;
  std::map<PVertex, IVertex> m_boundary_ivertices;
  const Planar_shape_type m_merge_type;
  const bool m_verbose;

  /*******************************
  **        MERGE PFACES        **
  ********************************/

  void merge_coplanar_pfaces(
    const std::size_t support_plane_idx) {

    CGAL_assertion_msg(false, "TODO: DELETE THIS ONE!");
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

    CGAL_assertion_msg(false, "TODO: DELETE THIS ONE!");
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

    CGAL_assertion_msg(false, "TODO: DELETE THIS ONE!");
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
  bool is_pface_inside_bbox(
    const std::size_t support_plane_idx,
    const std::vector<Point_2>& merged) const {

    CGAL_assertion_msg(false, "TODO: DELETE THIS ONE!");
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

    CGAL_assertion_msg(false, "TODO: DELETE THIS ONE!");
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

    std::set<IVertex> ivertices;
    const FT ptol = KSR::point_tolerance<FT>();
    for (const auto& iedge : iedges) {
      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);
      CGAL_assertion(isource != itarget);
      CGAL_assertion(KSR::distance(
        m_data.point_3(isource), m_data.point_3(itarget) ) >= ptol);
      ivertices.insert(isource);
      ivertices.insert(itarget);
    }

    for (const auto& ivertex : ivertices) {
      const auto point = m_data.to_2d(sp_idx, ivertex);
      points.push_back(std::make_pair(point,
        std::make_pair(m_data.null_pvertex(), ivertex)));
    }

    CGAL_assertion(points.size() == (pvertices.size() + ivertices.size()));
    // std::cout << "- num unique 1: " << points.size() << std::endl;

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
    //   // std::cout <<
    //   // m_data.str(pair.second.first) << " : " <<
    //   // m_data.str(pair.second.second) << std::endl;
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
      vh->info().sp_idx = pface.first;

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
        // std::cout <<
        // m_data.str(data.first) << " : " <<
        // m_data.str(data.second) << std::endl;
        // std::cout << m_data.to_3d(sp_idx, point) << std::endl;
        CGAL_assertion(idx < pedge_map.size());
        pedge_map[idx].push_back(data.second);
      }
    }

    CGAL_assertion(pedge_map.size() == n);
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& psource = polygon[i];
      const auto& ptarget = polygon[ip];

      CGAL_assertion(psource != ptarget);
      CGAL_assertion(vhs_pv.find(psource) != vhs_pv.end());
      CGAL_assertion(vhs_pv.find(ptarget) != vhs_pv.end());

      const auto vh_source = vhs_pv.at(psource);
      const auto vh_target = vhs_pv.at(ptarget);
      CGAL_assertion(vh_source != vh_target);
      CGAL_assertion(KSR::distance(
        vh_source->point(), vh_target->point()) >= ptol);

      m_input.insert(psource);
      if (pedge_map[i].size() > 0) {
        if (pedge_map[i].size() == 1 && sp_idx >= 6) {
          // CGAL_assertion_msg(false, "TODO: ADD TWO INTERMEDIATE IEDGES!");
          // In this case, we are sure, we do not have iedges along polygon edges,
          // so we simply insert missing constraints.

          const auto& iv = pedge_map[i][0];
          CGAL_assertion(vhs_iv.find(iv) != vhs_iv.end());
          const auto vh_mid = vhs_iv.at(iv);
          m_cdt.insert_constraint(vh_source, vh_mid);
          m_cdt.insert_constraint(vh_mid, vh_target);
          vh_mid->info().is_boundary_ivertex = true;
          // print_edge("original1", pface.first, vh_source, vh_mid);
          // print_edge("original2", pface.first, vh_mid, vh_target);

        } else if (pedge_map[i].size() > 1 && sp_idx >= 6) {
          CGAL_assertion_msg(false, "TODO: ADD MULTIPLE INTERMEDIATE IEDGES!");
          // If this case ever happens, we need to find out, which iedges from
          // all iedges (inserted as constraints below), are actually inserted
          // and if necessary add missing constraints connecting these iedges
          // to the polygon vertices, e.g. /pv/--add--/iv/--iedge--/iv/--edge--/pv/.
        } else {
          CGAL_assertion_msg(sp_idx < 6, "ERROR: WRONG CONSTRAINT CASE!");
        }
      } else {
        // CGAL_assertion_msg(false, "TODO: ADD STANDARD CONSTRAINT CASE!");
        // In this case, we do not have any intermediate ivertices along the pedge,
        // so we simply insert this pedge as a constraint.
        m_cdt.insert_constraint(vh_source, vh_target);
        // print_edge("original", pface.first, vh_source, vh_target);
      }
    }

    // Set pedge indices.
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t im = (i + n - 1) % n;
      const auto& pvertex = polygon[i];
      CGAL_assertion(vhs_pv.find(pvertex) != vhs_pv.end());
      const auto vh = vhs_pv.at(pvertex);
      CGAL_assertion(vh->info().pedge_indices.size() == 0);
      vh->info().pedge_indices.push_back(im);
      vh->info().pedge_indices.push_back(i);
    }

    for (std::size_t i = 0; i < pedge_map.size(); ++i) {
      const auto& pedge_ivertices = pedge_map[i];
      if (pedge_ivertices.size() == 0) continue;
      for (const auto& ivertex : pedge_ivertices) {
        CGAL_assertion(vhs_iv.find(ivertex) != vhs_iv.end());
        const auto vh = vhs_iv.at(ivertex);
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
    // std::cout << "num iedges: " << iedges.size() << std::endl;
    for (const auto& iedge : iedges) {
      const auto isource = m_data.source(iedge);
      const auto itarget = m_data.target(iedge);

      CGAL_assertion(isource != itarget);
      CGAL_assertion(vhs_iv.find(isource) != vhs_iv.end());
      CGAL_assertion(vhs_iv.find(itarget) != vhs_iv.end());

      const auto vh_source = vhs_iv.at(isource);
      const auto vh_target = vhs_iv.at(itarget);
      CGAL_assertion(vh_source != vh_target);
      CGAL_assertion(KSR::distance(
        vh_source->point(), vh_target->point()) >= ptol);

      const auto cid = m_cdt.insert_constraint(vh_source, vh_target);
      CGAL_assertion(m_map_intersections.find(cid) == m_map_intersections.end());
      m_map_intersections.insert(std::make_pair(cid, iedge));
    }

    // std::cout << "- num cdt verts 3: " << m_cdt.number_of_vertices() << std::endl;
    // std::cout << "- num cdt faces 3: " << m_cdt.number_of_faces()    << std::endl;

    // Add all points, which are not in unique points but in cdt.
    for (auto vit = m_cdt.finite_vertices_begin();
    vit != m_cdt.finite_vertices_end(); ++vit) {

      if (vit->info().sp_idx != KSR::no_element()) continue;
      if (vit->info().pvertex != m_data.null_pvertex()) continue;
      if (vit->info().ivertex != m_data.null_ivertex()) continue;

      const auto& point = vit->point();
      CGAL_assertion(vit->info().pvertex == m_data.null_pvertex());
      CGAL_assertion(vit->info().ivertex == m_data.null_ivertex());
      CGAL_assertion(!is_pvertex(vhs_pv, polygon, point));
      CGAL_assertion(!is_ivertex(pface.first, iedges, point));

      // std::cout << m_data.to_3d(sp_idx, point) << std::endl;
      const std::size_t idx = find_pedge(vhs_pv, polygon, point);
      // std::cout << "found idx: " << idx << std::endl;
      if (idx != KSR::no_element()) {
        CGAL_assertion(idx < polygon.size());
        CGAL_assertion(vit->info().pedge_indices.size() == 0);
        vit->info().pedge_indices.push_back(idx);
      }
      vit->info().sp_idx = pface.first;
    }
  }

  bool is_pvertex(
    const std::map<PVertex, Vertex_handle>& vhs_pv,
    const std::vector<PVertex>& polygon,
    const Point_2& query) const {

    const FT ptol = KSR::point_tolerance<FT>();
    for (const auto& pvertex : polygon) {
      CGAL_assertion(vhs_pv.find(pvertex) != vhs_pv.end());
      const auto vh = vhs_pv.at(pvertex);
      const auto& point = vh->point();
      const FT distance = KSR::distance(point, query);
      if (distance < ptol) return true;
    }
    return false;
  }

  bool is_ivertex(
    const std::size_t sp_idx,
    const std::set<IEdge>& iedges,
    const Point_2& query) {

    std::set<IVertex> ivertices;
    for (const auto& iedge : iedges) {
      ivertices.insert(m_data.source(iedge));
      ivertices.insert(m_data.target(iedge));
    }
    CGAL_assertion(ivertices.size() > 0);

    const FT ptol = KSR::point_tolerance<FT>();
    for (const auto& ivertex : ivertices) {
      const auto point = m_data.to_2d(sp_idx, ivertex);
      const FT distance = KSR::distance(point, query);
      if (distance < ptol) return true;
    }
    return false;
  }

  std::size_t find_pedge(
    const std::map<PVertex, Vertex_handle>& vhs_pv,
    const std::vector<PVertex>& polygon,
    const Point_2& query) const {

    const FT tol = KSR::tolerance<FT>();
    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& psource = polygon[i];
      const auto& ptarget = polygon[ip];

      CGAL_assertion(psource != ptarget);
      CGAL_assertion(vhs_pv.find(psource) != vhs_pv.end());
      CGAL_assertion(vhs_pv.find(ptarget) != vhs_pv.end());

      const auto vh_source = vhs_pv.at(psource);
      const auto vh_target = vhs_pv.at(ptarget);
      CGAL_assertion(vh_source != vh_target);

      const auto& source = vh_source->point();
      const auto& target = vh_target->point();
      CGAL_assertion(KSR::distance(source, target) >= KSR::point_tolerance<FT>());

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

    CGAL_assertion(face_index > 0);
    // std::cout << "- number of interior pfaces: " << face_index << std::endl;
  }

  bool is_boundary(const Edge& edge) const {

    const auto& fh = edge.first;
    const std::size_t idx = edge.second;
    const auto vh1 = fh->vertex( (idx + 1) % 3 );
    const auto vh2 = fh->vertex( (idx + 2) % 3 );

    if (m_cdt.is_infinite(vh1) || m_cdt.is_infinite(vh2)) {
      return false;
    }
    CGAL_assertion(!m_cdt.is_infinite(vh1));
    CGAL_assertion(!m_cdt.is_infinite(vh2));

    if (!m_cdt.is_constrained(edge)) {
      // print_edge("f0", vh1, vh2);
      return false;
    }

    const auto& pes1 = vh1->info().pedge_indices;
    const auto& pes2 = vh2->info().pedge_indices;
    CGAL_assertion(pes1.size() <= 2);
    CGAL_assertion(pes2.size() <= 2);

    if (pes1.size() == 0 || pes2.size() == 0) {
      // print_edge("f1", vh1, vh2);
      return false;
    }
    CGAL_assertion(pes1.size() > 0);
    CGAL_assertion(pes2.size() > 0);

    for (const std::size_t pe1 : pes1) {
      for (const std::size_t pe2 : pes2) {
        if (pe1 == pe2) {
          // print_edge("t0", vh1, vh2);
          return true;
        }
      }
    }
    // print_edge("f2", vh1, vh2);
    return false;
  }

  void print_edge(
    const std::string name, const std::size_t sp_idx,
    const Vertex_handle vh1, const Vertex_handle vh2) const {

    CGAL_assertion(sp_idx != KSR::no_element());
    std::cout << name << ": ";
    std::cout << m_data.to_3d(sp_idx, vh1->point()) << " ";
    std::cout << m_data.to_3d(sp_idx, vh2->point()) << std::endl;
  }

  void print_edge(
    const std::string name, const Vertex_handle vh1, const Vertex_handle vh2) const {

    CGAL_assertion(vh1->info().sp_idx != KSR::no_element());
    CGAL_assertion(vh2->info().sp_idx != KSR::no_element());
    std::cout << name << ": ";
    std::cout << m_data.to_3d(vh1->info().sp_idx, vh1->point()) << " ";
    std::cout << m_data.to_3d(vh2->info().sp_idx, vh2->point()) << std::endl;
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

          // Handle ivertices on the polygon boundary.
          if (source->info().is_boundary_ivertex) {
            CGAL_assertion(source->info().ivertex != m_data.null_ivertex());
            m_boundary_ivertices[source->info().pvertex] = source->info().ivertex;
          }
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

    CGAL_assertion(num_pfaces > 0);
    if (m_verbose) {
      std::cout << "** number of newly inserted pfaces: " << num_pfaces << std::endl;
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
        }
      }

      // Several incident intersections.
      // These are intersections of several iedges.
      if (is_frozen) {

        // Boundary ivertices.
        // These are not frozen since they are on the polygon boundary.
        if (m_boundary_ivertices.size() > 0) {
          const auto pit = m_boundary_ivertices.find(pvertex);
          if (pit != m_boundary_ivertices.end()) {
            const auto& pair = *pit;
            const auto& ivertex = pair.second;
            CGAL_assertion(pvertex == pair.first);
            set_boundary_ivertex(pvertex, ivertex);
            continue;
          }
        }

        // Interior ivertices. Frozen pvertex.
        CGAL_assertion(m_data.has_ivertex(pvertex));
        m_data.direction(pvertex) = CGAL::NULL_VECTOR;
        continue;
      }

      // No incident intersections = keep initial direction.
      // These are polygon vertices.
      if (iedge == m_data.null_iedge()) {
        CGAL_assertion(m_data.direction(pvertex) != CGAL::NULL_VECTOR);
        continue;
      }

      // Set future direction.
      // These are newly inserted points along the polygon boundary, which are
      // intersection points of multiple constraints. The simply follow the given iedge.
      m_data.connect(pvertex, iedge);
      const auto neighbors = get_polygon_neighbors(pvertex);
      Point_2 future_point; Vector_2 future_direction;
      compute_future_point_and_direction(
        pvertex, IVertex(), iedge,
        neighbors.first, neighbors.second, future_point, future_direction);
      CGAL_assertion(future_direction != Vector_2());
      m_data.direction(pvertex) = future_direction;
    }
  }

  void set_boundary_ivertex(
    const PVertex& pvertex, const IVertex& ivertex) {

    if (m_verbose) {
      std::cout.precision(20);
      std::cout << "*** setting boundary ivertex " << m_data.str(ivertex) <<
      " via pvertex " << m_data.str(pvertex) << std::endl;
      std::cout << "- ivertex: " << m_data.point_3(ivertex) << std::endl;
      std::cout << "- pvertex: " << m_data.point_3(pvertex) << std::endl;
    }

    // Get prev and next pvertices.
    PVertex prev, next;
    std::tie(prev, next) = m_data.border_prev_and_next(pvertex);
    if (m_verbose) {
      std::cout << "- prev: " << m_data.point_3(prev) << std::endl;
      std::cout << "- next: " << m_data.point_3(next) << std::endl;
    }

    // Freeze pvertex.
    const std::size_t sp_idx = pvertex.first;
    CGAL_assertion(sp_idx != KSR::no_element());
    m_data.direction(pvertex) = CGAL::NULL_VECTOR;
    const Point_2 ipoint = m_data.point_2(sp_idx, ivertex);
    m_data.support_plane(sp_idx).set_point(pvertex.second, ipoint);
    m_data.connect(pvertex, ivertex);

    // Get all connected iedges.
    std::vector< std::pair<IEdge, Direction_2> > iedges;
    m_data.get_and_sort_all_connected_iedges(sp_idx, ivertex, iedges);

    // Set reference directions.
    const auto prev_p = m_data.point_2(prev);
    const auto next_p = m_data.point_2(next);
    const Direction_2 ref_direction_prev(prev_p - ipoint);
    const Direction_2 ref_direction_next(next_p - ipoint);

    // Find the first iedge.
    std::size_t first_idx = std::size_t(-1);
    const std::size_t n = iedges.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;

      const auto& i_dir  = iedges[i].second;
      const auto& ip_dir = iedges[ip].second;
      if (ref_direction_next.counterclockwise_in_between(i_dir, ip_dir)) {
        first_idx = ip; break;
      }
    }
    CGAL_assertion(first_idx != std::size_t(-1));
    // std::cout << "- curr: " << m_data.segment_3(iedges[first_idx].first) << std::endl;

    // Find all crossed iedges.
    std::vector< std::pair<IEdge, bool> > crossed_iedges;
    std::size_t iedge_idx = first_idx;
    std::size_t iteration = 0;
    while (true) {
      const auto& iedge = iedges[iedge_idx].first;
      // std::cout << "- next: " << m_data.segment_3(iedge) << std::endl;

      if (iteration == iedges.size()) {
        CGAL_assertion_msg(iedges.size() == 2,
        "ERROR: SET BOUNDARY IVERTEX, CAN WE HAVE THIS CASE?");
        break;
      }

      const auto& ref_direction = iedges[iedge_idx].second;
      if (!ref_direction.counterclockwise_in_between(
        ref_direction_next, ref_direction_prev)) {
        break;
      }

      crossed_iedges.push_back(std::make_pair(iedge, false));
      iedge_idx = (iedge_idx + 1) % n;
      if (iteration >= iedges.size()) {
        CGAL_assertion_msg(false,
        "ERROR: SET BOUNDARY IVERTEX, WHY SO MANY ITERATIONS?");
      } ++iteration;
    }

    CGAL_assertion(crossed_iedges.size() >= 2);
    if (m_verbose) {
      std::cout << "- crossed " << crossed_iedges.size() << " iedges: " << std::endl;
      for (const auto& crossed_iedge : crossed_iedges) {
        std::cout << m_data.str(crossed_iedge.first) << ": " <<
        m_data.segment_3(crossed_iedge.first) << std::endl;
      }
    }

    // Compute future points and directions.
    std::vector<Point_2> future_points(2);
    std::vector<Vector_2> future_directions(2);
    const auto neighbors = get_polygon_neighbors(pvertex);
    compute_future_point_and_direction(
      pvertex, ivertex, crossed_iedges.front().first, neighbors.first, neighbors.second,
      future_points.front(), future_directions.front());
    compute_future_point_and_direction(
      pvertex, ivertex, crossed_iedges.back().first, neighbors.first, neighbors.second,
      future_points.back(), future_directions.back());

    // Crop the pvertex.
    std::vector<PVertex> new_pvertices;
    new_pvertices.resize(crossed_iedges.size(), m_data.null_pvertex());

    { // first crop
      if (m_verbose) std::cout << "- first crop" << std::endl;
      const auto cropped = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, next.second));
      CGAL_assertion(cropped != m_data.null_pvertex());

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      CGAL_assertion(cropped != pvertex);
      new_pvertices.front() = cropped;

      m_data.connect(pedge, crossed_iedges.front().first);
      m_data.connect(cropped, crossed_iedges.front().first);

      CGAL_assertion(future_directions.front() != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_points.front());
      m_data.direction(cropped) = future_directions.front();
      if (m_verbose) std::cout << "- cropped 1: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
    }

    { // second crop
      if (m_verbose) std::cout << "- second crop" << std::endl;
      const auto cropped = PVertex(pvertex.first, m_data.support_plane(pvertex).split_edge(pvertex.second, prev.second));
      CGAL_assertion(cropped != m_data.null_pvertex());

      const PEdge pedge(pvertex.first, m_data.support_plane(pvertex).edge(pvertex.second, cropped.second));
      CGAL_assertion(cropped != pvertex);
      new_pvertices.back() = cropped;

      m_data.connect(pedge, crossed_iedges.back().first);
      m_data.connect(cropped, crossed_iedges.back().first);

      CGAL_assertion(future_directions.back() != Vector_2());
      m_data.support_plane(cropped).set_point(cropped.second, future_points.back());
      m_data.direction(cropped) = future_directions.back();
      if (m_verbose) std::cout << "- cropped 2: " <<
        m_data.str(cropped) << ", " << m_data.point_3(cropped) << std::endl;
    }

    // Create new pfaces if any.
    m_data.add_pfaces(
      pvertex, ivertex, neighbors.first, neighbors.second,
      true, false, false, crossed_iedges, new_pvertices);

    // CGAL_assertion_msg(false, "TODO: HANDLE BOUNDARY IVERTICES!");
  }

  const std::pair<PVertex, PVertex> get_polygon_neighbors(
    const PVertex& pvertex) const {

    PVertex n1 = m_data.null_pvertex();
    PVertex n2 = m_data.null_pvertex();
    std::tie(n1, n2) = get_neighbors(pvertex);

    bool is_n1_okay = false;
    if (n1 != m_data.null_pvertex()) {
      is_n1_okay = update_neighbor(pvertex, n1);
    }

    bool is_n2_okay = false;
    if (n2 != m_data.null_pvertex()) {
      is_n2_okay = update_neighbor(pvertex, n2);
    }

    if (is_n1_okay && is_n2_okay) { } else {
      CGAL_assertion(is_n1_okay && !is_n2_okay);
      n2 = pvertex;
    }

    return std::make_pair(n1, n2);
  }

  const std::pair<PVertex, PVertex> get_neighbors(const PVertex& pvertex) const {

    std::pair<PVertex, PVertex> neighbors(
      m_data.null_pvertex(), m_data.null_pvertex());
    const auto pedges = m_data.pedges_around_pvertex(pvertex);
    for (const auto pedge : pedges) {
      // std::cout << "pedge: 2 " << m_data.segment_3(pedge) << " : "
      // << m_data.has_iedge(pedge) << std::endl;

      if (!m_data.has_iedge(pedge)) {
        const auto opposite = m_data.opposite(pedge, pvertex);

        if (neighbors.first == m_data.null_pvertex()) {
          neighbors.first = opposite;
          // std::cout << "assigned first neighbor: " << m_data.point_3(opposite) << std::endl;
        } else {
          CGAL_assertion(neighbors.first  != m_data.null_pvertex());
          CGAL_assertion(neighbors.second == m_data.null_pvertex());
          neighbors.second = opposite;
          // std::cout << "assigned second neighbor: " << m_data.point_3(opposite) << std::endl;
          break;
        }
      }
    }
    return neighbors;
  }

  // Set neighbor to the closest polygon vertex with the well-defined direction.
  bool update_neighbor(
    const PVertex& pvertex, PVertex& neighbor) const {

    bool is_found = (m_input.find(neighbor) != m_input.end());
    auto last = pvertex;
    auto curr = neighbor;
    while (!is_found) {
      PVertex next, ignored;
      std::tie(next, ignored) = m_data.border_prev_and_next(curr);
      if (next == last) {
        std::swap(next, ignored);
      }
      CGAL_assertion(ignored == last);

      last = curr; curr = next;
      if (m_input.find(curr) != m_input.end()) {
        neighbor = curr;
        is_found = true;
      }
    }
    return is_found;
  }

  void compute_future_point_and_direction(
    const PVertex& pvertex,
    const IVertex& ivertex, const IEdge& iedge,
    const PVertex& n1, const PVertex& n2,
    Point_2& future_point, Vector_2& future_direction) const {

    const std::size_t sp_idx = pvertex.first;
    CGAL_assertion(sp_idx != KSR::no_element());
    CGAL_assertion(sp_idx == n1.first);
    CGAL_assertion(sp_idx == n2.first);

    CGAL_assertion_msg(KSR::distance(
      m_data.point_2(sp_idx, m_data.source(iedge)),
      m_data.point_2(sp_idx, m_data.target(iedge))) >= KSR::point_tolerance<FT>(),
    "TODO: SET FUTURE DIRECTION, HANDLE ZERO-LENGTH IEDGE!");

    const bool is_debug = false;
    m_data.set_verbose(is_debug);
    const auto is_parallel = m_data.compute_future_point_and_direction(
      pvertex, ivertex, n1, n2, iedge, future_point, future_direction);
    m_data.set_verbose(m_verbose);

    CGAL_assertion_msg(!is_parallel,
    "TODO: COMPUTE FUTURE POINT AND DIRECTION, ADD PARALLEL CASE!");
    CGAL_assertion(future_direction != Vector_2());

    // std::cout << "curr point: " << m_data.point_3(pvertex) << std::endl;
    // std::cout << "futr point: " << m_data.to_3d(pvertex.first, future_point) << std::endl;
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
