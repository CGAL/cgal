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
// Author(s): Shepard Liu       <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_HALFEDGE_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_HALFEDGE_H

#include <algorithm>
#include <array>
#include <cstdlib>
#include <optional>
#include <type_traits>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/enum.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_has.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/** @brief Functor to approximate an x-monotone curve within an bounding box.
 *
 * The Approximation is done from xmin to xmax with a given step. For parts outbound the y limits and precedes or
 * succeeds a part within, the approximation may be skipped but there will be at least one point outside the bbox
 * for indication.
 */
template <typename Arrangement>
class Arr_bounded_approximate_halfedge {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Gt_point = typename Geom_traits::Point_2;

  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point = typename Approx_traits::Point;
  using Polyline = typename Approx_traits::Polyline;
  using Approx_kernel = typename Approx_traits::Approx_kernel;
  using Approx_line_2 = typename Approx_kernel::Line_2;
  using X_monotone_curve_2 = typename Geom_traits::X_monotone_curve_2;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;
  using Boundary_lines = std::array<Approx_line_2, 4>;

  constexpr static bool Has_approximate_xcv_with_bounds =
    has_approximate_xcv_with_bounds_v<Geom_traits, typename Geom_traits::Approximate_2>;

private:
  struct Context : public Bounded_render_context {
    Context(const Bounded_render_context& ctx, const X_monotone_curve_2& curve, Polyline& polyline) :
      Bounded_render_context(ctx),
      m_polyline(polyline), m_curve(curve)
    {}

    // Prevent accidental copying.
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

  public:
    /*! \brief Insert a point to the polyline if it is within the x-range of the curve
     * \note Will be replaced after AosApproximateUnboundedTraits_2 is fully available.
     * \param pt
     */
    void insert(Point pt) {
      if (pt.x() < this->xmin()) {
        // We need the last point if not yet x-inbound.
        m_last_pt = pt;
        return;
      }
      else if (pt.x() > this->xmax()) return;

      m_polyline.push_back(pt);
      m_last_pt = pt;
    }

    const std::optional<Point>& last_pt() const { return m_last_pt; }

  private:
    std::optional<Point> m_last_pt;

  public:
    Polyline& m_polyline;
    const X_monotone_curve_2& m_curve;
  };

  /*! \brief Computes the intersection point between the given boundary side and the line segment from last_pt to pt.
   */
  Point boundary_intersection(const Context& ctx, Point pt, Boundary_side side) const {
    std::optional<double> x, y;
    const Approx_line_2* line = nullptr;
    switch(side) {
    case Boundary_side::Left:
      x = ctx.xmin();
      line = &m_left;
      break;
    case Boundary_side::Right:
      x = ctx.xmax();
      line = &m_right;
      break;
    case Boundary_side::Top:
      y = ctx.ymax();
      line = &m_top;
      break;
    case Boundary_side::Bottom:
      y = ctx.ymin();
      line = &m_bottom;
      break;
    default:
      CGAL_assertion(false && "Unexpected side of boundary.");
    }
    Point inter = std::get<Point>(*CGAL::intersection(Approx_line_2(*ctx.last_pt(), pt), *line));
    if (x.has_value()) return Point(*x, inter.y());
    return Point(inter.x(), *y);
  }

  /*! \brief Trace approximated curve point in ltr ordering, adding boundary intersections if necessary.
   *
   * \note This method will eventually be replaced by AosApproximateUnboundedTraits_2.
   */
  void trace_add(Context& ctx, Point pt) const {
    if (! ctx.last_pt().has_value()) {
      ctx.insert(pt);
      return;
    }
    if (ctx.last_pt()->x() < ctx.xmin() && pt.x() >= ctx.xmin())
      ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Left));
    if (ctx.last_pt()->y() < ctx.ymin()) {
      if (pt.y() > ctx.ymin()) ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Bottom));
      if (pt.y() > ctx.ymax()) ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Top));
    }
    else if (ctx.last_pt()->y() > ctx.ymax()) {
      if (pt.y() < ctx.ymax()) ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Top));
      if (pt.y() < ctx.ymin()) ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Bottom));
    }
    else {
      if (pt.y() < ctx.ymin())
        ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Bottom));
      else if (pt.y() > ctx.ymax())
        ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Top));
    }
    if (ctx.last_pt()->x() <= ctx.xmax() && pt.x() > ctx.xmax())
      ctx.insert(boundary_intersection(ctx, pt, Boundary_side::Right));
    ctx.insert(pt);
  }

  /*! \brief Check if the point is within the x-range of the curve.
   */
  static bool is_in_x_range(const Context& ctx, const Gt_point& pt) {
    const Geom_traits& traits = ctx.m_traits;
    const X_monotone_curve_2& curve = ctx.m_curve;

    if constexpr(has_is_in_x_range_v<Geom_traits>) return curve.is_in_x_range(pt);
    if constexpr(!has_parameter_space_in_x_2<Geom_traits>::value) {
      const auto& min_pt = traits.construct_min_vertex_2_object()(curve);
      const auto& max_pt = traits.construct_max_vertex_2_object()(curve);
      return ((traits.compare_x_2_object()(pt, min_pt) != CGAL::SMALLER) &&
              (traits.compare_x_2_object()(pt, max_pt) != CGAL::LARGER));
    }

    Comparison_result left_cmp;
    if (auto left_loc = traits.parameter_space_in_x_2_object()(curve, ARR_MIN_END); left_loc == ARR_INTERIOR)
      left_cmp = traits.compare_x_2_object()(pt, traits.construct_min_vertex_2_object()(curve));
    else if (left_loc == ARR_LEFT_BOUNDARY)
      left_cmp = CGAL::LARGER;
    else
      left_cmp = traits.compare_x_on_boundary_2_object()(pt, curve, ARR_MIN_END);
    if (left_cmp == CGAL::SMALLER) return false;
    if (left_cmp == CGAL::EQUAL) return true;

    Comparison_result right_cmp;
    if (auto right_loc = traits.parameter_space_in_x_2_object()(curve, ARR_MAX_END); right_loc == ARR_INTERIOR)
      right_cmp = traits.compare_x_2_object()(pt, traits.construct_max_vertex_2_object()(curve));
    else if (right_loc == ARR_RIGHT_BOUNDARY)
      right_cmp = CGAL::SMALLER;
    else
      right_cmp = traits.compare_x_on_boundary_2_object()(pt, curve, ARR_MAX_END);
    return right_cmp != CGAL::LARGER;
  }

  /*! \brief transform approximated curve points(ltr ordering) in place based on the halfedge, giving correct
   * ordering, continuity, etc.
   */
  static void transform_polyline(Context& ctx, Polyline& polyline, const Halfedge_const_handle& he)
  { transform_polyline_impl<Geom_traits>(ctx, polyline, he); }

  // For planar arrangements, we only need to reverse the polyline if the halfedge is rtl.
  template <typename Gt, std::enable_if_t<!is_or_derived_from_curved_surf_traits_v<Gt>, int> = 0>
  static void transform_polyline_impl(Context&, Polyline& polyline, const Halfedge_const_handle& he) {
    if (he->direction() == CGAL::ARR_LEFT_TO_RIGHT) return;
    std::reverse(polyline.begin(), polyline.end());
  }

  template <typename Gt, std::enable_if_t<is_or_derived_from_agas_v<Gt>, int> = 0>
  static void transform_polyline_impl(Context& ctx, Polyline& polyline, const Halfedge_const_handle& he) {
    using Direction_3 = typename Geom_traits::Direction_3;
    using Vector_3 = typename Geom_traits::Vector_3;

    if (polyline.size() < 2) return;
    const X_monotone_curve_2& curve = he->curve();
    const auto& traits = ctx.m_traits;
    if (curve.is_vertical()) {
      Direction_3 normal_dir = curve.is_directed_right() ? curve.normal() : -curve.normal();
      Direction_3 azimuth_dir(CGAL::cross_product(Vector_3(0, 0, 1), normal_dir.vector()));
      Approx_nt azimuth = ctx.to_uv(traits.approximate_2_object()(traits.construct_point_2_object()(azimuth_dir))).x();
      if (azimuth == 0 && he->direction() == ARR_LEFT_TO_RIGHT) azimuth = 2 * CGAL_PI;
      std::transform(polyline.begin(), polyline.end(), polyline.begin(),
                     [azimuth](Point pt) { return Point(azimuth, pt.y()); });
    }
    else if (polyline.back().x() == 0) {
      // For strictly x-monotone arcs whose target point sits on the boundary, the x should be set to 2 * CGAL_PI
      polyline.back() = Point(2 * CGAL_PI, polyline.back().y());
    }
    if (he->direction() == CGAL::ARR_LEFT_TO_RIGHT) return;
    std::reverse(polyline.begin(), polyline.end());
  }

  void approximate_curve(Context& ctx) const { approximate_curve_impl<Geom_traits>(ctx); }

  // If Approximate_2 supports curve approximation with bounding box
  template <typename Gt, std::enable_if_t<has_approximate_xcv_with_bounds_v<Gt, typename Gt::Approximate_2>, int> = 0>
  void approximate_curve_impl(Context& ctx) const {
    const Geom_traits& traits = ctx.m_traits;
    const X_monotone_curve_2& curve = ctx.m_curve;
    Polyline& polyline = ctx.m_polyline;
    auto compare_y_at_x_2 = traits.compare_y_at_x_2_object();

    if (is_in_x_range(ctx, m_top_left)) {
      if (compare_y_at_x_2(m_top_left, curve) == CGAL::SMALLER) {
        polyline.insert(polyline.end(), {Approx_traits::Null_point, Point(ctx.xmin(), ctx.ymax())});
      }
      else if (compare_y_at_x_2(m_bottom_left, curve) == CGAL::LARGER) {
        polyline.insert(polyline.end(), {Approx_traits::Null_point, Point(ctx.xmin(), ctx.ymin())});
      }
    }
    traits.approximate_2_object()(curve, ctx.m_approx_error,
                                  boost::make_function_output_iterator([&ctx, this](Approx_point approx_pt)
                                  { ctx.m_polyline.push_back(snap_to_boundary(ctx, ctx.to_uv(approx_pt))); }),
                                  ctx.bbox(), true);
    if (is_in_x_range(ctx, m_top_right)) {
      if (compare_y_at_x_2(m_top_right, curve) == CGAL::SMALLER) {
        polyline.insert(polyline.end(), {Point(ctx.xmax(), ctx.ymax()), Approx_traits::Null_point});
      }
      else if (compare_y_at_x_2(m_bottom_right, curve) == CGAL::LARGER) {
        polyline.insert(polyline.end(), {Point(ctx.xmax(), ctx.ymin()), Approx_traits::Null_point});
      }
    }
  }

  // If Approximate_2 does not support curve approximation with bounding box
  template <typename Gt, std::enable_if_t<!has_approximate_xcv_with_bounds_v<Gt, typename Gt::Approximate_2>, int> = 0>
  void approximate_curve_impl(Context& ctx) const {
    auto approx = m_ctx.m_traits.approximate_2_object();
    approx(ctx.m_curve, ctx.m_approx_error,
           boost::make_function_output_iterator([&ctx, this](Approx_point pt) { trace_add(ctx, ctx.to_uv(pt)); }), true);
  }

  /*! \brief Adjusts a point by snapping it to the nearest boundary to reduce floating-point error.
   *
   * \return The adjusted (snapped) point if it lies within snapping tolerance, or the original point otherwise.
   */
  Point snap_to_boundary(const Context& ctx, Point pt) const {
    Approx_nt x = pt.x(), y = pt.y();
    if (std::abs(x - ctx.xmin()) < m_ep_left) x = ctx.xmin();
    else if (std::abs(x - ctx.xmax()) < m_ep_right) x = ctx.xmax();
    if (std::abs(y - ctx.ymin()) < m_ep_bottom) y = ctx.ymin();
    else if (std::abs(y - ctx.ymax()) < m_ep_top) y = ctx.ymax();
    return Point(x, y);
  }

public:
  Arr_bounded_approximate_halfedge(const Bounded_render_context& ctx) :
    m_ctx(ctx),
    m_left(ctx.bottom_left(), ctx.top_left()),
    m_right(ctx.bottom_right(), ctx.top_right()),
    m_bottom(ctx.bottom_left(), ctx.bottom_right()),
    m_top(ctx.top_left(), ctx.top_right()) {
    Construct_gt_point_2<Geom_traits> ctr_p;
    m_top_left = ctr_p(ctx.to_cartesian(ctx.top_left()));
    m_top_right = ctr_p(ctx.to_cartesian(ctx.top_right()));
    m_bottom_left = ctr_p(ctx.to_cartesian(ctx.bottom_left()));
    m_bottom_right = ctr_p(ctx.to_cartesian(ctx.bottom_right()));
    Approx_nt ep_base = std::numeric_limits<Approx_nt>::epsilon();
    m_ep_left = std::max(std::abs(ep_base * ctx.xmin()), ep_base);
    m_ep_right = std::max(std::abs(ep_base * ctx.xmax()), ep_base);
    m_ep_bottom = std::max(std::abs(ep_base * ctx.ymin()), ep_base);
    m_ep_top = std::max(std::abs(ep_base * ctx.ymax()), ep_base);
  }

  const Polyline& operator()(const Halfedge_const_handle& he) const {
    CGAL_assertion(!he->is_fictitious());

    auto& cache = m_ctx.m_cache.halfedges();
    auto [iter, inserted] = cache.try_emplace(he, Polyline());
    Polyline& polyline = iter->second;
    if (!inserted) return polyline;
    if (m_ctx.is_cancelled()) return polyline;

    const X_monotone_curve_2& curve = he->curve();
    Context ctx(m_ctx, curve, polyline);
    approximate_curve(ctx);
    Polyline poly_copy(polyline);
    transform_polyline(ctx, polyline, he);

    // also approximate the twin halfedge
    auto [twin_iter, twin_inserted] = cache.try_emplace(he->twin(), std::move(poly_copy));
    if (twin_inserted) transform_polyline(ctx, twin_iter->second, he->twin());
    // The previous iterator might have been invalidated by the second try_emplace call, so we do an extra lookup.
    return cache.at(he);
  }

private:
  const Bounded_render_context& m_ctx;
  Approx_line_2 m_left, m_right, m_bottom, m_top;
  Gt_point m_top_left, m_top_right, m_bottom_left, m_bottom_right;
  Approx_nt m_ep_left, m_ep_right, m_ep_bottom, m_ep_top;
};

} // namespace draw_aos
} // namespace CGAL

#endif
