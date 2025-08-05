#ifndef CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H
#define CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
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
#include <CGAL/draw_constrained_triangulation_2.h>

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
  using Point_geom = typename Approx_traits::Point_geom;
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
  // For planar arrangements, Constrained_triangulation_2 is enough.
  using Ct = std::conditional_t<is_or_derived_from_curved_surf_traits<Geom_traits>,
                                Constrained_Delaunay_triangulation_2<Epick, Tds, Exact_predicates_tag>,
                                Constrained_triangulation_2<Epick, Tds, Exact_predicates_tag>>;
  using KPoint = Epick::Point_2;
  using KPoint_with_info = std::pair<KPoint, Point_index>;

  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

public:
  using Insert_iterator = boost::function_output_iterator<std::function<void(Point_geom)>>;

private:
  static KPoint transform_point(Point_geom pt) { return KPoint(pt.x(), pt.y()); }

  static Point_geom offset_boundary_point(Point_geom pt, Side_of_boundary side, double offset) {
    CGAL_precondition(side != Side_of_boundary::None);
    switch(side) {
    case Side_of_boundary::Left:
      return Point_geom(pt.x() - offset, pt.y());
    case Side_of_boundary::Right:
      return Point_geom(pt.x() + offset, pt.y());
    case Side_of_boundary::Top:
      return Point_geom(pt.x(), pt.y() + offset);
    case Side_of_boundary::Bottom:
      return Point_geom(pt.x(), pt.y() - offset);
    default:
      return pt; // Should not reach here
    }
  }

  /**
   * Sample the interior of a face for geometry traits on curved surfaces
   */
  void sample_interior() { sample_interior_impl<Geom_traits>(); }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  void sample_interior_impl() {}

  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  void sample_interior_impl() {
    const double scale_factor = 20.0;
    double cell_size = m_ctx.m_approx_error * scale_factor;
    std::size_t grid_x_min = m_face_bbox.xmin() / cell_size;
    std::size_t grid_x_max = std::ceil(m_face_bbox.xmax() / cell_size);
    std::size_t grid_y_min = m_face_bbox.ymin() / cell_size;
    std::size_t grid_y_max = std::ceil(m_face_bbox.ymax() / cell_size);
    for(std::size_t x = grid_x_min + 1; x < grid_x_max; ++x)
      for(std::size_t y = grid_y_min + 1; y < grid_y_max; ++y)
        m_points.emplace_back(Point_geom(x * cell_size, y * cell_size));
    // sample the bbox boundary
    for(double x = m_face_bbox.xmin(); x <= m_face_bbox.xmax(); x += cell_size) {
      m_points.emplace_back(Point_geom(x, m_face_bbox.ymin()));
      m_points.emplace_back(Point_geom(x, m_face_bbox.ymax()));
    }
    for(double y = m_face_bbox.ymin(); y <= m_face_bbox.ymax(); y += cell_size) {
      m_points.emplace_back(Point_geom(m_face_bbox.xmin(), y));
      m_points.emplace_back(Point_geom(m_face_bbox.xmax(), y));
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
    // Constrained_triangulation_2 and Constrained_Delaunay_triangulation_2 have slightly different interfaces.
    if constexpr(is_or_derived_from_curved_surf_traits<Geom_traits>)
      m_ct.insert(transformed_begin, transformed_end);
    else
      m_ct.template insert_with_info<KPoint_with_info>(transformed_begin, transformed_end);
  }

  Side_of_boundary shared_boundary(const Point_geom& pt1, const Point_geom& pt2) const {
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

  void add_boundary_helper_point(Point_geom from, Point_geom to) {
    if(from == to) return;
    auto shared_side = shared_boundary(from, to);
    if(shared_side == Side_of_boundary::None) return;
    m_helper_indices.push_back(m_points.size());
    m_points.emplace_back(
        offset_boundary_point(Point_geom((from.x() + to.x()) / 2, (from.y() + to.y()) / 2), shared_side, m_offset));
    m_offset += 0.5;
  }

public:
  Arr_bounded_face_triangulator(const Bounded_render_context& ctx)
      : m_ctx(ctx) {}

  Insert_iterator insert_iterator() {
    return boost::make_function_output_iterator(std::function([this](Point_geom pt) {
      CGAL_assertion_msg(m_ctx.contains(pt), "Point should be within the bounding box.");

      if(m_points.size() != 0) add_boundary_helper_point(m_points.back(), pt);
      m_points.emplace_back(pt);
      this->m_face_bbox += pt.bbox();
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

    m_constraint_end = m_points.size();
    add_boundary_helper_point(m_points.back(), m_points.front());
    sample_interior();
    insert_ccb();
    // insert_constraint() should be called after insert_with_info(), or info will not be set correctly.
    m_ct.insert_constraint(boost::make_transform_iterator(m_points.begin(), transform_point),
                           boost::make_transform_iterator(m_points.begin() + m_constraint_end, transform_point), true);
    if(m_ct.number_of_faces() == 0) {
      return Triangulated_face();
    }

#if defined(CGAL_DRAW_AOS_DEBUG)
    debug_print(*this);
#endif

    unordered_flat_map<typename Ct::Face_handle, bool> in_domain_map;
    in_domain_map.reserve(m_ct.number_of_faces());
    boost::associative_property_map<decltype(in_domain_map)> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(m_ct, in_domain);
    // Collect triangles within the constrained domain.
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
  std::vector<Point_geom> m_points;
  Bbox_2 m_face_bbox;                        // The bounding box of the face.
  std::vector<std::size_t> m_helper_indices; // The extra helper point indices when inserting outer ccb constraint.
  std::size_t m_constraint_end;              // The index past the last point in the outer ccb constraint.
  double m_offset{0.5};
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
  for(std::size_t i = 0; i < triangulator.m_constraint_end; ++i) {
    const auto& pt = m_points[i];
    if(helper_indices_index < m_helper_indices.size() && i == m_helper_indices[helper_indices_index]) {
      helper_indices_index++;
      continue;
    }
    ofs_ccb << pt.x() << " " << pt.y() << "\n";
  }

  std::ofstream ofs_ccb_constraint(ccb_constraint_file_path);
  for(std::size_t i = 0; i < triangulator.m_constraint_end; ++i) {
    const auto& pt = m_points[i];
    ofs_ccb_constraint << pt.x() << " " << pt.y() << "\n";
  }
}
#endif

} // namespace draw_aos
} // namespace CGAL

#endif