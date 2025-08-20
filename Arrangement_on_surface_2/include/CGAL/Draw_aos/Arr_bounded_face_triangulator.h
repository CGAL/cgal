// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shepard Liu	 <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H
#define CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H

#include <algorithm>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/iterator/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>

#if defined(CGAL_DRAW_AOS_DEBUG) && defined(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR)
#include <fstream>
#include <filesystem>

template <typename Arrangement>
class Arr_bounded_face_triangulator;

template <typename Arrangement>
void debug_print(const Arr_bounded_face_triangulator<Arrangement>& triangulator);
#endif

namespace CGAL {
namespace draw_aos {

/**
 * @brief Triangulator for a face of an arrangement within a bounding box.
 */
template <typename Arrangement>
class Arr_bounded_face_triangulator
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  constexpr static bool Is_on_curved_surface = is_or_derived_from_curved_surf_traits_v<Geom_traits>;

  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Point = typename Approx_traits::Point;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approx_kernel = typename Approx_traits::Approx_kernel;
  using Triangle_soup = typename Approx_traits::Triangle_soup;
  using Triangle = typename Triangle_soup::Triangle;
  using Face_const_handle = typename Arrangement::Face_const_handle;

#if defined(CGAL_DRAW_AOS_DEBUG)
  template <typename T>
  friend void debug_print(const Arr_bounded_face_triangulator<T>& triangulator);
#endif

  enum Point_type { Vertex_only, Constraint_only, Vertex_and_constraint };

  /*!
   * \brief A index wrapper defaulted to invalid.
   */
  class Index
  {
  public:
    Index() = default;
    Index(int idx)
        : m_index(idx) {}

    bool is_valid() const { return m_index != Invalid_index; }
    operator int() const { return m_index; }

  private:
    constexpr static int Invalid_index = -1;
    int m_index{Invalid_index};
  };

  using Epick = Exact_predicates_inexact_constructions_kernel;
  using Vb = Triangulation_vertex_base_with_info_2<Index, Epick>;
  using Fb = Constrained_triangulation_face_base_2<Epick>;
  using Tds = Triangulation_data_structure_2<Vb, Fb>;
  // For planar arrangements, Constrained_triangulation_2 is enough.
  using Ct = std::conditional_t<Is_on_curved_surface,
                                Constrained_Delaunay_triangulation_2<Epick, Tds, Exact_predicates_tag>,
                                Constrained_triangulation_2<Epick, Tds, Exact_predicates_tag>>;

  using KPoint = Epick::Point_2;
  using KPoint_with_index = std::pair<KPoint, Index>;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

public:
  using value_type = Point;

private:
  static KPoint to_kpoint(Point pt) { return KPoint(pt.x(), pt.y()); }

  /*!
   * \brief Offset a point on a specific boundary outward by a given offset.
   *
   * \pre side != Boundary_side::None
   */
  static Point offset_boundary_point(Point pt, Boundary_side side, double offset) {
    CGAL_precondition(side != Boundary_side::None);

    switch(side) {
    case Boundary_side::Left:
      return Point(pt.x() - offset, pt.y());
    case Boundary_side::Right:
      return Point(pt.x() + offset, pt.y());
    case Boundary_side::Top:
      return Point(pt.x(), pt.y() + offset);
    case Boundary_side::Bottom:
      return Point(pt.x(), pt.y() - offset);
    default:
      return pt; // Should not reach here
    }
  }

  /*!
   * \brief Find the shared boundary side of two points, or None if they are not on the same boundary.
   */
  Boundary_side shared_boundary(const Point& pt1, const Point& pt2) const {
    if(m_ctx.is_on_left(pt1) && m_ctx.is_on_left(pt2)) return Boundary_side::Left;
    if(m_ctx.is_on_right(pt1) && m_ctx.is_on_right(pt2)) return Boundary_side::Right;
    if(m_ctx.is_on_bottom(pt1) && m_ctx.is_on_bottom(pt2)) return Boundary_side::Bottom;
    if(m_ctx.is_on_top(pt1) && m_ctx.is_on_top(pt2)) return Boundary_side::Top;
    return Boundary_side::None;
  }

  /*!
   * \brief Add a helper point on the shared boundary of two points if they are on the same boundary side.
   *
   * When triangulating a arrangement face within a bounding box, curves outside the bounding box are projected on the
   * four sides of the bbox. Topological errors could be introduced if several polylines are lying on the same side.
   * Thus we add the midpoint in between the two points on boundary and move it outward with an increasing offset.
   */
  void add_boundary_helper_point(Point from, Point to) {
    // Arrangements on curved surfaces currently draws the entire parameter space, so there's no need to add
    // helper points.
    if constexpr(Is_on_curved_surface) return;
    if(from == to) return;
    auto shared_side = shared_boundary(from, to);
    if(shared_side == Boundary_side::None) return;
    Point mid = CGAL::midpoint(from, to);
    m_points.push_back(offset_boundary_point(mid, shared_side, m_offset += 0.1));
    m_point_types.push_back(Constraint_only);
  }

  void insert_all_vertices() {
    auto vertex_filter = [this](int idx) { return m_point_types[idx] != Constraint_only; };
    auto index_to_point_with_info = [this](int idx) -> KPoint_with_index {
      return std::make_pair(to_kpoint(m_points[idx]), idx);
    };
    auto indexes_begin = boost::make_counting_iterator<int>(0);
    auto indexes_end = boost::make_counting_iterator<int>(m_points.size());
    auto filtered_begin = boost::make_filter_iterator(vertex_filter, indexes_begin, indexes_end);
    auto filtered_end = boost::make_filter_iterator(vertex_filter, indexes_end, indexes_end);
    auto transformed_begin = boost::make_transform_iterator(filtered_begin, index_to_point_with_info);
    auto transformed_end = boost::make_transform_iterator(filtered_end, index_to_point_with_info);

    // Constrained_triangulation_2 and Constrained_Delaunay_triangulation_2 have slightly different interfaces.
    if constexpr(Is_on_curved_surface)
      m_ct.insert(transformed_begin, transformed_end);
    else
      m_ct.template insert_with_info<KPoint_with_index>(transformed_begin, transformed_end);
  }

  void insert_all_constraints() {
    auto constraint_filter = [this](int idx) { return m_point_types[idx] != Vertex_only; };
    auto index_to_point = [this](int idx) -> KPoint { return to_kpoint(m_points[idx]); };
    for(auto [start_idx, end_idx] : m_cst_ranges) {
      auto indexes_begin = boost::make_counting_iterator<int>(start_idx);
      auto indexes_end = boost::make_counting_iterator<int>(end_idx);
      auto filtered_begin = boost::make_filter_iterator(constraint_filter, indexes_begin, indexes_end);
      auto filtered_end = boost::make_filter_iterator(constraint_filter, indexes_end, indexes_end);
      auto transformed_begin = boost::make_transform_iterator(filtered_begin, index_to_point);
      auto transformed_end = boost::make_transform_iterator(filtered_end, index_to_point);
      m_ct.insert_constraint(transformed_begin, transformed_end, true);
    }
  }

public:
  Arr_bounded_face_triangulator(const Bounded_render_context& ctx, Face_const_handle fh)
      : m_ctx(ctx)
      , m_fh(fh) {}

  void push_back(Point pt) {
    CGAL_assertion_msg(m_curr_cst_begin.has_value(), "Call start_constraint() before push_back().");

    if(m_points.size() - *m_curr_cst_begin >= 2) add_boundary_helper_point(m_points.back(), pt);
    m_points.push_back(pt);
    m_point_types.push_back(Vertex_and_constraint);
  }

  void start_constraint() { m_curr_cst_begin = m_points.size(); }

  void end_constraint() {
    CGAL_assertion_msg(m_curr_cst_begin.has_value(), "Call start_constraint() before end_constraint().");

    int cst_begin = *m_curr_cst_begin;
    m_curr_cst_begin.reset();
    if(m_points.size() - cst_begin <= 2) {
      m_points.erase(m_points.begin() + cst_begin, m_points.end());
      m_point_types.erase(m_point_types.begin() + cst_begin, m_point_types.end());
      return;
    }
    add_boundary_helper_point(m_points.back(), m_points[cst_begin]);
    m_cst_ranges.emplace_back(cst_begin, m_points.size());
  }

  /*!
   * \brief Converts the triangulator to a triangulated face, moving internal data to the result.
   *
   * \return Triangulated_face
   */
  operator Triangle_soup() && {
    CGAL_assertion_msg(!m_curr_cst_begin.has_value(), "Call end_constraint() before conversion");

    if(m_points.empty()) return Triangle_soup();
    if constexpr(Is_on_curved_surface) {
      if(auto it = m_ctx.m_face_points.find(m_fh); it != m_ctx.m_face_points.end()) {
        m_points.insert(m_points.end(), it->second.begin(), it->second.end());
        m_point_types.insert(m_point_types.end(), it->second.size(), Vertex_only);
      }
    }
    insert_all_vertices();
    insert_all_constraints();
    if(m_ct.number_of_faces() == 0) return Triangle_soup();

#if defined(CGAL_DRAW_AOS_DEBUG)
    debug_print(*this);
#endif

    unordered_flat_map<typename Ct::Face_handle, bool> in_domain_map;
    in_domain_map.reserve(m_ct.number_of_faces());
    boost::associative_property_map<decltype(in_domain_map)> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(m_ct, in_domain);
    // Collect triangles within the constrained domain.
    Triangle_soup ts;
    ts.triangles.reserve(m_ct.number_of_faces());
    for(auto fit = m_ct.finite_faces_begin(); fit != m_ct.finite_faces_end(); ++fit) {
      Index v1 = fit->vertex(0)->info();
      Index v2 = fit->vertex(1)->info();
      Index v3 = fit->vertex(2)->info();
      if(!v1.is_valid() || !v2.is_valid() || !v3.is_valid()) continue;
      if(!get(in_domain, fit)) continue;
      ts.triangles.push_back(Triangle{v1, v2, v3});
    }
    ts.points = std::move(m_points);
    return ts;
  }

private:
  const Bounded_render_context& m_ctx;
  Face_const_handle m_fh;
  Ct m_ct;
  std::vector<Point> m_points;
  std::vector<Point_type> m_point_types;
  std::vector<std::pair<int, int>> m_cst_ranges;
  std::optional<int> m_curr_cst_begin;
  double m_offset{0};
};

#if defined(CGAL_DRAW_AOS_DEBUG) && defined(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR)
template <typename Arrangement>
void debug_print(const Arr_bounded_face_triangulator<Arrangement>& triangulator) {
  const auto& ctx = triangulator.m_ctx;
  const auto& m_points = triangulator.m_points;
  const auto& m_point_types = triangulator.m_point_types;
  using Point_type = typename Arr_bounded_face_triangulator<Arrangement>::Point_type;

  using Path = std::filesystem::path;
  Path debug_dir(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR);
  std::string index_file_name = "index.txt";
  Path index_file_path = debug_dir / index_file_name;
  std::string points_file_name_prefix = "face_" + std::to_string(*ctx.debug_counter) + "_points";
  std::string ccb_constraint_file_name_prefix = "face_" + std::to_string(*ctx.debug_counter) + "_constraint";
  const_cast<int&>(*ctx.debug_counter)++;

  std::ofstream ofs_index(index_file_path, std::ios::app);

  auto points_filename = points_file_name_prefix + ".txt";
  auto points_path = debug_dir / points_filename;
  std::ofstream ofs_points(points_path);
  ofs_index << points_filename << std::endl;
  for(int i = 0; i < triangulator.m_points.size(); ++i) {
    if(m_point_types[i] == Point_type::Constraint_only) continue;
    const auto& pt = m_points[i];
    ofs_points << pt.x() << " " << pt.y() << "\n";
  }

  int counter = 0;
  for(auto [start_idx, end_idx] : triangulator.m_cst_ranges) {
    auto filename = ccb_constraint_file_name_prefix + "_" + std::to_string(counter++) + ".txt";
    auto filepath = debug_dir / filename;
    ofs_index << filename << std::endl;
    std::ofstream ofs_ccb_constraint(filepath);
    for(int i = start_idx; i < end_idx; ++i) {
      if(m_point_types[i] == Point_type::Vertex_only) continue;
      const auto& pt = m_points[i];
      ofs_ccb_constraint << pt.x() << " " << pt.y() << "\n";
    }
  }
}
#endif

} // namespace draw_aos
} // namespace CGAL

#endif