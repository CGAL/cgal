#ifndef CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H
#define CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H

#include "CGAL/Constrained_triangulation_2.h"
#include "CGAL/Constrained_triangulation_face_base_2.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include "CGAL/Draw_aos/Arr_render_context.h"
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Triangulation_vertex_base_with_info_2.h"
#include "CGAL/mark_domain_in_triangulation.h"
#include "CGAL/number_utils.h"
#include "CGAL/unordered_flat_map.h"
#include <cstddef>
#include <CGAL/Draw_aos/helpers.h>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <utility>
#include <vector>

namespace CGAL {

/**
 * @brief Triangulator a bounded face of an arrangement
 *
 * @note The face must have an outer CCB.
 */
class Arr_bounded_face_triangulator
{
  using Approx_geom_traits = Arr_approximation_geometry_traits;
  using Approx_kernel = Approx_geom_traits::Approximation_kernel;
  using Approx_point = Approx_geom_traits::Approx_point;
  using Point_vec = Approx_geom_traits::Apporx_point_vec;
  using Triangle = Approx_geom_traits::Triangle;
  using Triangulated_face = Approx_geom_traits::Triangulated_face;

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
  using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Vb = CGAL::Triangulation_vertex_base_with_info_2<Point_index, Epick>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<Epick>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  using Ct = Constrained_triangulation_2<Epick, Tds, Exact_predicates_tag>;
  using KPoint = Epick::Point_2;
  using KPoint_with_info = std::pair<KPoint, Point_index>;

  /**
   * @brief RAII-style inserter for one CCB in a triangulation.
   * Collects points and inserts them as a constraint on destruction.
   * Only one instance per Arr_face_triangulator is allowed at a time.
   */
  template <typename Ccb_tag>
  class Ccb_constraint
  {
    constexpr static bool Is_outer_ccb = std::is_same_v<Ccb_tag, Outer_ccb_tag>;

    friend class Arr_bounded_face_triangulator;
    using Side_of_boundary = Arr_bounded_render_context::Side_of_boundary;

    Ccb_constraint(Arr_bounded_face_triangulator& triangulator)
        : m_triangulator(&triangulator)
        , m_ccb_start(m_triangulator->m_points.size()) {
      triangulator.m_has_active_constraint = true;
    }

  private:
    std::vector<KPoint_with_info>& points() { return m_triangulator->m_points; }
    const std::vector<KPoint_with_info>& points() const { return m_triangulator->m_points; }
    auto first_point() const { return points()[m_ccb_start].first; }
    auto last_point() const { return points().back().first; }
    const Arr_bounded_render_context& ctx() const { return m_triangulator->m_ctx; }
    void add_point(const KPoint& pt) { points().emplace_back(pt, points().size()); }
    std::size_t ccb_size() const { return m_triangulator->m_points.size() - m_ccb_start - m_helper_indices.size(); }

    KPoint offset_boundary_point(const KPoint& pt, Side_of_boundary side, double offset) const {
      CGAL_precondition(side != Side_of_boundary::None);
      switch(side) {
      case Side_of_boundary::Left:
        return KPoint(pt.x() - offset, pt.y());
      case Side_of_boundary::Right:
        return KPoint(pt.x() + offset, pt.y());
      case Side_of_boundary::Top:
        return KPoint(pt.x(), pt.y() + offset);
      case Side_of_boundary::Bottom:
        return KPoint(pt.x(), pt.y() - offset);
      default:
        return pt; // Should not reach here
      }
    }

    void insert_ccb() {
      auto& ct = m_triangulator->m_ct;
      auto& points = m_triangulator->m_points;

      auto begin = points.begin() + m_ccb_start;
      auto end = points.end();

      if constexpr(Is_outer_ccb) {
        auto helpers_iter = m_helper_indices.begin();
        auto helpers_end = m_helper_indices.end();

        auto concrete_pt_filter = [&helpers_iter, helpers_end](std::size_t idx) {
          if(helpers_iter == helpers_end) {
            return true;
          }
          if(idx == *helpers_iter) {
            ++helpers_iter;
            return false;
          }
          return true;
        };

        auto index_to_point_with_info = [&points](std::size_t idx) { return points[idx]; };
        auto indexes_begin = boost::make_counting_iterator<std::size_t>(m_ccb_start);
        auto indexes_end = boost::make_counting_iterator<std::size_t>(points.size());
        auto filtered_begin = boost::make_filter_iterator(concrete_pt_filter, indexes_begin, indexes_end);
        auto filtered_end = boost::make_filter_iterator(concrete_pt_filter, indexes_end, indexes_end);
        auto transformed_begin = boost::make_transform_iterator(filtered_begin, index_to_point_with_info);
        auto transformed_end = boost::make_transform_iterator(filtered_end, index_to_point_with_info);

        ct.insert_with_info<KPoint_with_info>(transformed_begin, transformed_end);
      } else {
        ct.insert_with_info<KPoint_with_info>(begin, end);
      }
    }

    void try_add_offset(const KPoint& from, const KPoint& to) {
      if(from == to) {
        return;
      }
      auto shared_side = ctx().shared_boundary_side(from, to);
      if(shared_side == Arr_bounded_render_context::Side_of_boundary::None) {
        return;
      }
      // m_helper_indices.push_back(points().size());
      // add_point(offset_boundary_point(from, shared_side, m_offset));
      m_helper_indices.push_back(points().size());
      add_point(offset_boundary_point(KPoint((from.x() + to.x()) / 2, (from.y() + to.y()) / 2), shared_side, m_offset));

      // TODO: we'll come back to find out if we do need to offset the second point.
      // m_helper_indices.push_back(points().size());
      // add_point(offset_boundary_point(to, shared_side, m_offset));
      // Doesn't matter how much we increase the offset.
      m_offset += 0.5;
    }

  public:
    Ccb_constraint(const Ccb_constraint&) = delete;
    Ccb_constraint& operator=(const Ccb_constraint&) = delete;
    Ccb_constraint(Ccb_constraint&& other) noexcept
        : m_triangulator(other.m_triangulator) {
      other.m_triangulator = nullptr;
    }
    Ccb_constraint& operator=(Ccb_constraint&& other) noexcept {
      if(this != &other) {
        m_triangulator = other.m_triangulator;
        other.m_triangulator = nullptr;
      }
      return *this;
    }

    auto insert_iterator() {
      return boost::make_function_output_iterator([this](const Approx_point& pt) {
        CGAL_assertion_msg(m_triangulator != nullptr, "Use of destructed or moved Ccb_constraint object.");
        CGAL_assertion_msg(ctx().contains(pt), "Outbound point in Ccb_constraint.");
        KPoint kp(pt.x(), pt.y());

        if(Is_outer_ccb && ccb_size() != 0) {
          try_add_offset(last_point(), kp);
        }

        add_point(kp);
      });
    }

    ~Ccb_constraint() {
      if(m_triangulator == nullptr) {
        return;
      }

      if(Is_outer_ccb && ccb_size() != 0) {
        try_add_offset(last_point(), first_point());
      }

      insert_ccb();
      m_triangulator->m_has_active_constraint = false;
      m_triangulator = nullptr;
    }

  private:
    Arr_bounded_face_triangulator* m_triangulator;
    const std::size_t m_ccb_start;

    // These are used only for outer CCBs.
    double m_offset = 0.5;                     // Doesn't matter how much we offset.
    std::vector<std::size_t> m_helper_indices; // The offseted point indices when inserting outer ccb constraint
  };

public:
  Arr_bounded_face_triangulator(const Arr_bounded_render_context& ctx)
      : m_ctx(ctx) {}
  Arr_bounded_face_triangulator(const Arr_bounded_face_triangulator& other) = delete;
  Arr_bounded_face_triangulator& operator=(const Arr_bounded_face_triangulator& other) = delete;

  operator Triangulated_face() && {
    CGAL_assertion(!m_has_active_constraint && "There is an active constraint in the triangulator.");
    CGAL_assertion(m_outer_ccb_processed && "Outer CCB has not been processed yet.");
    if(m_points.empty() || m_ct.number_of_faces() == 0) {
      return Triangulated_face();
    }

    // insert_constraint() should be called after insert_with_info(), or info will not be set correctly.
    auto first_of_pair = [](const KPoint_with_info& pt_with_info) { return pt_with_info.first; };
    for(std::size_t i = 0; i < m_ccb_start_indices.size(); ++i) {
      auto begin = m_points.begin() + m_ccb_start_indices[i];
      auto end = i + 1 < m_ccb_start_indices.size() ? m_points.begin() + m_ccb_start_indices[i + 1] : m_points.end();

      m_ct.insert_constraint(boost::make_transform_iterator(begin, first_of_pair),
                             boost::make_transform_iterator(end, first_of_pair), true);
    }

    unordered_flat_map<Ct::Face_handle, bool> in_domain_map;
    boost::associative_property_map<decltype(in_domain_map)> in_domain(in_domain_map);
    CGAL::mark_domain_in_triangulation(m_ct, in_domain);

    Triangulated_face tf;
    tf.triangles.reserve(m_ct.number_of_faces());
    for(auto fit = m_ct.finite_faces_begin(); fit != m_ct.finite_faces_end(); ++fit) {
      if(!get(in_domain, fit)) {
        continue;
      }
      Point_index v1 = fit->vertex(0)->info();
      Point_index v2 = fit->vertex(1)->info();
      Point_index v3 = fit->vertex(2)->info();
      if(!v1.is_valid() || !v2.is_valid() || !v3.is_valid()) {
        continue;
      }
      Triangle tri{v1, v2, v3};
      tf.triangles.emplace_back(tri);
    }

    auto transform_first_of_pair = [](const KPoint_with_info& pt_with_info) {
      return Approx_point(CGAL::to_double(pt_with_info.first.x()), CGAL::to_double(pt_with_info.first.y()));
    };
    tf.points = Point_vec(boost::make_transform_iterator(m_points.begin(), transform_first_of_pair),
                          boost::make_transform_iterator(m_points.end(), transform_first_of_pair));
    return tf;
  }

  /**
   * @brief Get an constraint object for a certain ccb.
   * @note Only one Ccb_constraint can be active per triangulator at a certain time point. When it goes out of scope,
   * the points collected will be inserted as a constraint into the triangulation.
   */
  template <typename Ccb_tag>
  Ccb_constraint<Ccb_tag> make_ccb_constraint() {
    if constexpr(std::is_same_v<Ccb_tag, Outer_ccb_tag>) {
      CGAL_assertion_msg(!m_outer_ccb_processed, "Outer CCB has already been processed.");
      m_outer_ccb_processed = true;
    }
    CGAL_assertion_msg(!m_has_active_constraint, "Only one Ccb_constraint can be active per triangulator at a time.");
    m_ccb_start_indices.push_back(m_points.size());
    return Ccb_constraint<Ccb_tag>(*this);
  }

private:
  const Arr_bounded_render_context& m_ctx;
  Ct m_ct;
  std::vector<KPoint_with_info> m_points;
  std::vector<std::size_t> m_ccb_start_indices;
  bool m_has_active_constraint = false;
  bool m_outer_ccb_processed = false;
};

} // namespace CGAL

#endif // CGAL_DRAW_AOS_ARR_FACE_TRIANGULATOR_H