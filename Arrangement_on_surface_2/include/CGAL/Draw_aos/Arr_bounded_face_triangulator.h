#ifndef CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H
#define CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H

#include <algorithm>
#include <cstddef>
#include <functional>
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
#include <CGAL/Draw_aos/Arr_render_context.h>

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
 *
 * The original topology of a face is preserved, but the geometry will be bounded by the bbox.
 */
template <typename Arrangement>
class Arr_bounded_face_triangulator
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point_vec = typename Approx_traits::Approx_point_vec;
  using Triangle = typename Approx_traits::Triangle;
  using Triangulated_face = typename Approx_traits::Triangulated_face;

#if defined(CGAL_DRAW_AOS_DEBUG)
  template <typename T>
  friend void debug_print(const Arr_bounded_face_triangulator<T>& triangulator);
#endif

  struct Point_index
  {
    constexpr static std::size_t Invalid_index = std::numeric_limits<std::size_t>::max();
    std::size_t index{Invalid_index};
    Point_index() = default;
    Point_index(std::size_t idx)
        : index(idx) {}
    bool is_valid() const { return index != Invalid_index; }
    operator std::size_t() const { return index; }
  };
  using Epick = Exact_predicates_inexact_constructions_kernel;
  using Vb = Triangulation_vertex_base_with_info_2<Point_index, Epick>;
  using Fb = Constrained_triangulation_face_base_2<Epick>;
  using Tds = Triangulation_data_structure_2<Vb, Fb>;
  using Ct = Constrained_triangulation_2<Epick, Tds, Exact_predicates_tag>;
  using KPoint = Epick::Point_2;
  using KPoint_with_info = std::pair<KPoint, Point_index>;

  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

public:
  using Insert_iterator = boost::function_output_iterator<std::function<void(Approx_point)>>;

private:
  static KPoint transform_point(Approx_point pt) { return KPoint(pt.x(), pt.y()); }

  static Approx_point offset_boundary_point(Approx_point pt, Side_of_boundary side, double offset) {
    CGAL_precondition(side != Side_of_boundary::None);
    switch(side) {
    case Side_of_boundary::Left:
      return Approx_point(pt.x() - offset, pt.y());
    case Side_of_boundary::Right:
      return Approx_point(pt.x() + offset, pt.y());
    case Side_of_boundary::Top:
      return Approx_point(pt.x(), pt.y() + offset);
    case Side_of_boundary::Bottom:
      return Approx_point(pt.x(), pt.y() - offset);
    default:
      return pt; // Should not reach here
    }
  }

  void insert_ccb() {
    auto begin = m_points.begin();
    auto end = m_points.end();
    auto helpers_iter = m_helper_indices.begin();
    auto helpers_end = m_helper_indices.end();
    // A two pointers algorithm implemented with boost filters.
    auto point_filter = [&helpers_iter, helpers_end](std::size_t idx) {
      if(helpers_iter == helpers_end) {
        return true;
      }
      if(idx == *helpers_iter) {
        ++helpers_iter;
        return false;
      }
      return true;
    };
    auto index_to_point_with_info = [this](std::size_t idx) -> KPoint_with_info {
      return std::make_pair(transform_point(m_points[idx]), idx);
    };
    auto indexes_begin = boost::make_counting_iterator<std::size_t>(0);
    auto indexes_end = boost::make_counting_iterator<std::size_t>(m_points.size());
    auto filtered_begin = boost::make_filter_iterator(point_filter, indexes_begin, indexes_end);
    auto filtered_end = boost::make_filter_iterator(point_filter, indexes_end, indexes_end);
    auto transformed_begin = boost::make_transform_iterator(filtered_begin, index_to_point_with_info);
    auto transformed_end = boost::make_transform_iterator(filtered_end, index_to_point_with_info);
    m_ct.template insert_with_info<KPoint_with_info>(transformed_begin, transformed_end);
  }

  Side_of_boundary shared_boundary(const Approx_point& pt1, const Approx_point& pt2) const {
    if(pt1.x() == m_ctx.xmin() && pt2.x() == m_ctx.xmin() && m_ctx.contains_y(pt1.y()) && m_ctx.contains_y(pt2.y()))
      return Side_of_boundary::Left;
    if(pt1.x() == m_ctx.xmax() && pt2.x() == m_ctx.xmax() && m_ctx.contains_y(pt1.y()) && m_ctx.contains_y(pt2.y()))
      return Side_of_boundary::Right;
    if(pt1.y() == m_ctx.ymin() && pt2.y() == m_ctx.ymin() && m_ctx.contains_x(pt1.x()) && m_ctx.contains_x(pt2.x()))
      return Side_of_boundary::Bottom;
    if(pt1.y() == m_ctx.ymax() && pt2.y() == m_ctx.ymax() && m_ctx.contains_x(pt1.x()) && m_ctx.contains_x(pt2.x()))
      return Side_of_boundary::Top;
    return Side_of_boundary::None;
  }

  void add_boundary_helper_point(Approx_point from, Approx_point to) {
    auto shared_side = shared_boundary(from, to);
    if(shared_side == Side_of_boundary::None) return;
    m_helper_indices.push_back(m_points.size());
    m_points.emplace_back(
        offset_boundary_point(Approx_point((from.x() + to.x()) / 2, (from.y() + to.y()) / 2), shared_side, m_offset));
    m_offset += 0.5; // It can be arbitrary increment.
  }

public:
  Arr_bounded_face_triangulator(const Bounded_render_context& ctx)
      : m_ctx(ctx) {}

  Insert_iterator insert_iterator() {
    return boost::make_function_output_iterator(std::function([this](Approx_point pt) {
      CGAL_assertion_msg(m_ctx.contains(pt), "Outbound point in Ccb_constraint.");

      if(m_points.size() != 0) {
        add_boundary_helper_point(m_points.back(), pt);
      }
      m_points.emplace_back(pt);
    }));
  }

  /**
   * @brief Converts the triangulator to a triangulated face, moving internal data to the result.
   *
   * @return Triangulated_face
   */
  operator Triangulated_face() && {
    if(m_points.empty()) {
      return Triangulated_face();
    }

    add_boundary_helper_point(m_points.back(), m_points.front());
    insert_ccb();
    // insert_constraint() should be called after insert_with_info(), or info will not be set correctly.
    m_ct.insert_constraint(boost::make_transform_iterator(m_points.begin(), transform_point),
                           boost::make_transform_iterator(m_points.end(), transform_point), true);

#if defined(CGAL_DRAW_AOS_DEBUG)
    debug_print(*this);
#endif

    if(m_ct.number_of_faces() == 0) {
      return Triangulated_face();
    }

    unordered_flat_map<typename Ct::Face_handle, bool> in_domain_map;
    in_domain_map.reserve(m_ct.number_of_faces());
    boost::associative_property_map<decltype(in_domain_map)> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(m_ct, in_domain);

    Triangulated_face tf;
    tf.triangles.reserve(m_ct.number_of_faces());
    for(auto fit = m_ct.finite_faces_begin(); fit != m_ct.finite_faces_end(); ++fit) {
      if(!get(in_domain, fit)) continue;
      Point_index v1 = fit->vertex(0)->info();
      Point_index v2 = fit->vertex(1)->info();
      Point_index v3 = fit->vertex(2)->info();
      if(!v1.is_valid() || !v2.is_valid() || !v3.is_valid()) continue;
      Triangle tri{v1, v2, v3};
      tf.triangles.emplace_back(tri);
    }
    tf.points = std::move(m_points);
    return tf;
  }

private:
  const Bounded_render_context& m_ctx;
  Ct m_ct;
  std::vector<Approx_point> m_points;
  double m_offset = 0.5;                     // Doesn't matter how much we offset.
  std::vector<std::size_t> m_helper_indices; // The offseted point indices when inserting outer ccb constraint
};

#if defined(CGAL_DRAW_AOS_DEBUG) && defined(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR)
template <typename Arrangement>
void debug_print(const Arr_bounded_face_triangulator<Arrangement>& triangulator) {
  const auto& ctx = triangulator.m_ctx;
  const auto& m_points = triangulator.m_points;
  const auto& m_helper_indices = triangulator.m_helper_indices;

  using Path = std::filesystem::path;
  Path debug_dir(CGAL_DRAW_AOS_TRIANGULATOR_DEBUG_FILE_DIR);
  std::string index_file_name = "index.txt";
  Path index_file_path = debug_dir / index_file_name;
  std::string ccb_file_name = "ccb_" + std::to_string(*ctx.debug_counter) + ".txt";
  std::string ccb_constraint_file_name = "ccb_constraint_" + std::to_string(*ctx.debug_counter) + ".txt";
  Path ccb_file_path = debug_dir / ccb_file_name;
  Path ccb_constraint_file_path = debug_dir / ccb_constraint_file_name;
  const_cast<int&>(*ctx.debug_counter)++;

  std::ofstream ofs_index(index_file_path, std::ios::app);
  ofs_index << ccb_file_name << "\n" << ccb_constraint_file_name << std::endl;

  std::ofstream ofs_ccb(ccb_file_path);
  std::size_t helper_indices_index = 0;
  for(std::size_t i = 0; i < m_points.size(); ++i) {
    const auto& pt = m_points[i];
    if(helper_indices_index < m_helper_indices.size() && i == m_helper_indices[helper_indices_index]) {
      helper_indices_index++;
      continue;
    }
    ofs_ccb << pt.x() << " " << pt.y() << "\n";
  }

  std::ofstream ofs_ccb_constraint(ccb_constraint_file_path);
  for(const auto& pt : m_points) {
    ofs_ccb_constraint << pt.x() << " " << pt.y() << "\n";
  }
}
#endif

} // namespace draw_aos
} // namespace CGAL

#endif