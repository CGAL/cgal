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

#include <CGAL/license/Arrangement_on_surface_2.h>

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
    /*! \brief Insert a point to the polyline if it is within the x-range of the bounding box.
     * A point to the left of the x-range is retained as the last point, so that the intersection with the left side
     * of the bounding box can be computed when the next point is traced.
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
   * It is used for both the approximation of a curve without constraints and the approximation of a curve clipped to
   * the bounding box (see `AosTraits::ApproximateUnbounded_2`). Points that lie below or above the bounding box are
   * kept; they are clamped onto the bottom or top side when faces are approximated.
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
    if (polyline.size() < 2) return;
    const X_monotone_curve_2& curve = he->curve();
    const auto& traits = ctx.m_traits;
    if (traits.is_vertical_2_object()(curve)) {
      // All points of a vertical curve share the same azimuth (x), except for
      // the poles, where the azimuth is undefined and the conversion sets it to
      // 0. Assign the same azimuth to all points.
      std::optional<Approx_nt> azimuth;
      if (traits.is_on_y_identification_2_object()(curve)) {
        // The curve lies on the identification curve, where the azimuth is
        // either 0 or 2 * pi depending on the side of the halfedge. Set it
        // exactly, as the approximated points may deviate slightly.
        azimuth = (he->direction() == ARR_LEFT_TO_RIGHT) ? 2 * CGAL_PI : 0;
      }
      else {
        // Obtain the azimuth from a point that is not a pole (0 < y < pi).
        auto it = std::find_if(polyline.begin(), polyline.end(),
                               [](Point pt) { return (0 < pt.y()) && (pt.y() < CGAL_PI); });
        CGAL_assertion_msg(it != polyline.end(),
                           "The approximation of a vertical curve has no point other than poles.");
        if (it != polyline.end()) azimuth = it->x();
      }
      if (azimuth.has_value())
        std::transform(polyline.begin(), polyline.end(), polyline.begin(),
                       [x = *azimuth](Point pt) { return Point(x, pt.y()); });
    }
    else if (polyline.back().x() == 0) {
      // For strictly x-monotone arcs whose target point sits on the boundary, the x should be set to 2 * CGAL_PI
      polyline.back() = Point(2 * CGAL_PI, polyline.back().y());
    }
    if (he->direction() == CGAL::ARR_LEFT_TO_RIGHT) return;
    std::reverse(polyline.begin(), polyline.end());
  }

  void approximate_curve(Context& ctx) const { approximate_curve_impl<Geom_traits>(ctx); }

  // If Approximate_2 supports curve approximation with bounding box.
  // The approximation is clipped to the bounding box, except that every
  // portion of the curve that lies below or above the box (within the
  // vertical strip of the box) is represented by points strictly below or
  // above the box. Clamping the points to the box, when faces are
  // approximated, turns these portions into the bottom or top side of the
  // box, which yields the corners of the box. Edges are drawn only between
  // points inside the box, so these portions are not drawn. The points are
  // traced exactly as the points of the approximation of a bounded curve.
  template <typename Gt, std::enable_if_t<has_approximate_xcv_with_bounds_v<Gt, typename Gt::Approximate_2>, int> = 0>
  void approximate_curve_impl(Context& ctx) const {
    auto approx = m_ctx.m_traits.approximate_2_object();
    approx(ctx.m_curve, ctx.m_approx_error,
           boost::make_function_output_iterator([&ctx, this](Approx_point pt) { trace_add(ctx, ctx.to_uv(pt)); }),
           ctx.bbox(), true);
  }

  // If Approximate_2 does not support curve approximation with bounding box
  template <typename Gt, std::enable_if_t<!has_approximate_xcv_with_bounds_v<Gt, typename Gt::Approximate_2>, int> = 0>
  void approximate_curve_impl(Context& ctx) const {
    auto approx = m_ctx.m_traits.approximate_2_object();
    approx(ctx.m_curve, ctx.m_approx_error,
           boost::make_function_output_iterator([&ctx, this](Approx_point pt) { trace_add(ctx, ctx.to_uv(pt)); }), true);
  }

public:
  Arr_bounded_approximate_halfedge(const Bounded_render_context& ctx) :
    m_ctx(ctx),
    m_left(ctx.bottom_left(), ctx.top_left()),
    m_right(ctx.bottom_right(), ctx.top_right()),
    m_bottom(ctx.bottom_left(), ctx.bottom_right()),
    m_top(ctx.top_left(), ctx.top_right())
  {}

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
};

} // namespace draw_aos
} // namespace CGAL

#endif
