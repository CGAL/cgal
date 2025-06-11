#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_CURVE_2_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_CURVE_2_H

#include "CGAL/Arr_enums.h"
#include "CGAL/Draw_aos/Arr_bounded_approximate_point_2.h"
#include "CGAL/Draw_aos/Arr_bounded_compute_y_at_x.h"
#include "CGAL/Draw_aos/Arr_construct_curve_end.h"
#include "CGAL/Draw_aos/Arr_construct_segments.h"
#include "CGAL/Draw_aos/Arr_render_context.h"
#include "CGAL/Draw_aos/helpers.h"
#include "CGAL/Draw_aos/Arr_approximation_geometry_traits.h"
#include "CGAL/basic.h"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <vector>

namespace CGAL {

/**
 * @brief Functor to approximate an x-monotone curve within an bounding box.
 * The bbox here has closed boundary.
 *
 * The Approximation is done from xmin to xmax with a given step. For parts outbound the y limits and precedes or
 * succeeds a part within, the approximation may be skipped but there will be at least one point outside the bbox for
 * indication.
 *
 * @note Bounded approximation is meaningful only when the curve has at least two points within the bbox (boundary
 * points included).
 *
 * TODO: Possible optimizations:
 * - Specialize for traits that models Approximate_2 on curves.
 */
class Arr_bounded_approximate_curve_2
{
  using FT = Geom_traits::FT;
  using Point_2 = Geom_traits::Point_2;
  using X_monotone_curve_2 = Geom_traits::X_monotone_curve_2;
  using Halfedge_const_handle = Arrangement::Halfedge_const_iterator;
  using Approx_point = Arr_approximation_geometry_traits::Approx_point;
  using Point_geom = Arr_approximation_geometry_traits::Point_geom;
  using Polyline_geom = Arr_approximation_geometry_traits::Polyline_geom;
  using Intersections_vector = std::vector<Point_2>;

  struct Execution_context : public Arr_context_delegator<Arr_bounded_render_context>
  {
    Execution_context(const Arr_bounded_render_context& ctx,
                      const X_monotone_curve_2& curve,
                      const Arr_bounded_approximate_point_2& approx_pt,
                      const Arr_bounded_compute_y_at_x& compute_y_at_x,
                      const Intersections_vector& top_inters,
                      const Intersections_vector& bottom_inters,
                      Polyline_geom& polyline)
        : Arr_context_delegator(ctx)
        , curve(curve)
        , bounded_compute_y_at_x(compute_y_at_x)
        , bounded_approx_pt(approx_pt)
        , top_inters(top_inters)
        , bottom_inters(bottom_inters)
        , min_end(Arr_construct_curve_end<Geom_traits>(ctx.traits)(curve, ARR_MIN_END))
        , max_end(Arr_construct_curve_end<Geom_traits>(ctx.traits)(curve, ARR_MAX_END))
        , tight_xmin(is_min_end_bounded() ? std::clamp(min_end->x(), FT(ctx.xmin()), FT(ctx.xmax())) : FT(ctx.xmin()))
        , tight_xmax(is_max_end_bounded() ? std::clamp(max_end->x(), FT(ctx.xmin()), FT(ctx.xmax())) : FT(ctx.xmax()))
        , tight_ymin(is_min_end_bounded() ? std::clamp(min_end->y(), FT(ctx.ymin()), FT(ctx.ymax())) : FT(ctx.ymin()))
        , tight_ymax(is_max_end_bounded() ? std::clamp(max_end->y(), FT(ctx.ymin()), FT(ctx.ymax())) : FT(ctx.ymax()))
        , out_it(std::back_inserter(polyline)) {}

    bool has_y_intersections() const { return !top_inters.empty() || !bottom_inters.empty(); }
    bool is_min_end_bounded() const { return min_end.has_value(); }
    bool is_max_end_bounded() const { return max_end.has_value(); }
    bool is_bounded_curve() const { return is_min_end_bounded() && is_max_end_bounded(); }

    const X_monotone_curve_2& curve;
    const Arr_bounded_compute_y_at_x& bounded_compute_y_at_x;
    const Arr_bounded_approximate_point_2& bounded_approx_pt;
    const Intersections_vector &top_inters, bottom_inters;
    const std::optional<Point_2> min_end, max_end;
    const FT tight_xmin, tight_xmax, tight_ymin, tight_ymax;
    std::back_insert_iterator<Polyline_geom> out_it;
  };

private:
  static std::vector<Point_2> compute_intersections(const X_monotone_curve_2& cv1,
                                                    const X_monotone_curve_2& cv2,
                                                    const typename Geom_traits::Intersect_2& intersect_2,
                                                    const Arr_construct_curve_end<Geom_traits>& cst_curve_end) {
    using Intersect_point = std::pair<Geom_traits::Point_2, Geom_traits::Multiplicity>;
    using Intersect_curve = Geom_traits::X_monotone_curve_2;
    using Intersect_type = std::variant<Intersect_point, Intersect_curve>;

    std::vector<Point_2> intersections;
    auto out_it = std::back_inserter(intersections);
    intersect_2(cv1, cv2, boost::make_function_output_iterator([&out_it, &cst_curve_end](const Intersect_type& res) {
                  if(auto* pt = std::get_if<Intersect_point>(&res)) {
                    *out_it++ = pt->first;
                    return;
                  }
                  if(auto* cv = std::get_if<Intersect_curve>(&res)) {
                    *out_it++ = cst_curve_end(*cv, ARR_MIN_END).value();
                    *out_it++ = cst_curve_end(*cv, ARR_MAX_END).value();
                    return;
                  }
                  CGAL_assertion(false && "Unexpected intersection type");
                }));
    return intersections;
  }

  static std::optional<Point_2> first_intersection(Execution_context& ctx) {
    if(!ctx.top_inters.empty() && !ctx.bottom_inters.empty()) {
      return ctx->compare_xy_2(ctx.top_inters.front(), ctx.bottom_inters.front()) == CGAL::SMALLER
                 ? ctx.top_inters.front()
                 : ctx.bottom_inters.front();
    }
    if(!ctx.top_inters.empty()) {
      return ctx.top_inters.front();
    }
    if(!ctx.bottom_inters.empty()) {
      return ctx.bottom_inters.front();
    }
    return std::nullopt;
  }

  /**
   * @brief approximate strictly x-monotone curve segment that does not cross the y bounds.
   *
   * @precondition: The segment is either inbound or outbound the bbox in the given range.
   * @param start the x-coordinate of the segment start(exclusive)
   * @param end the x-coordinate of the segment end(exclusive)
   * @param step the step to approximate the curve segment, negative values allowed.
   * @returns true if this part of the curve is within the closed bbox
   */
  static void approximate_simple_curve_segment(Execution_context& ctx, const FT& start, const FT& end, double step) {
    for(FT x = start + step; x < end; x += step) {
      auto y = ctx.bounded_compute_y_at_x(ctx.curve, x);
      if(!y.has_value()) {
        // break as soon as there's no more intersections
        break;
      }
      if(y == ctx->ymin() || y == ctx->ymax()) {
        // The segment overlaps with the bbox edge. There's no need to insert a dummy point.
        break;
      }
      *ctx.out_it++ = ctx->approx_pt(Point_2(x, y.value()));
      if(y > ctx->ymax() || y < ctx->ymin()) {
        // We are outside the bbox. The dummy point is inserted to indicate the curve is outside the bbox.
        break;
      }
    }
  };

  static void approximate_vertical_curve(Execution_context& ctx) {
    if(ctx.is_bounded_curve() && ctx.min_end->x() < ctx->xmin() && ctx.min_end->x() > ctx->xmax()) {
      // The curve is outside the bbox in x direction, no need to approximate
      return;
    }
    if(!ctx.is_bounded_curve() && !ctx.has_y_intersections()) {
      // The curve has unbounded end and has no intersections with the top or bottom edges,
      // it must be outbound in x direction.
      return;
    }
    // The vertical curve is now within the x bounds.

    FT tymin = ctx.tight_ymin;
    FT tymax = ctx.tight_ymax;
    if(tymax == tymin) {
      // But the curve is not degenerate. So, either it has only one point within the bbox or it is
      // entirely outside the bbox in y direction.
      return;
    }
    // Now we gaurantee that the curve has at least two points within the bbox in y direction.
    // We have to obtain the x coordinate of this vertical curve.
    FT x = ctx.is_bounded_curve() ? ctx.min_end->x() : first_intersection(ctx)->x();
    *ctx.out_it++ = ctx->approx_pt(Point_2(x, tymin));
    *ctx.out_it++ = ctx->approx_pt(Point_2(x, tymax));
  }

public:
  Arr_bounded_approximate_curve_2(const Arr_bounded_render_context& ctx,
                                  const Arr_bounded_approximate_point_2& point_approx)
      : m_bounded_compute_y_at_x(ctx)
      , m_approx_pt(point_approx)
      , m_ctx(ctx)
      , m_top(ctx.cst_horizontal_segment(ctx.ymax(), ctx.xmin(), ctx.xmax()))
      , m_bottom(ctx.cst_horizontal_segment(ctx.ymin(), ctx.xmin(), ctx.xmax())) {}

  /**
   * @brief Approximate an x-monotone curve from left to right within the bounding box.
   *
   * @param he non-fictitious halfedge handle
   * @return const Polyline_geom&
   */
  const Polyline_geom& operator()(const Halfedge_const_handle& he) const {
    CGAL_assertion(!he->is_fictitious());

    auto [polyline, inserted] = m_ctx.cache.try_emplace(he);
    if(!inserted) {
      return polyline;
    }

    if(m_ctx.is_cancelled()) {
      return polyline;
    }

    const X_monotone_curve_2& curve = he->curve();

    if(curve.is_degenerate()) {
      return polyline;
    }

    auto top_inters = compute_intersections(curve, m_top, m_ctx.intersect_2, m_ctx.cst_curve_end);
    auto bottom_inters = compute_intersections(curve, m_bottom, m_ctx.intersect_2, m_ctx.cst_curve_end);
    Execution_context ctx(m_ctx, curve, m_approx_pt, m_bounded_compute_y_at_x, top_inters, bottom_inters, polyline);

    if(ctx->is_vertical_2(curve)) {
      approximate_vertical_curve(ctx);
      return polyline;
    }

    polyline.reserve(top_inters.size() + bottom_inters.size());

    FT txmin = ctx.tight_xmin;
    FT txmax = ctx.tight_xmax;
    FT last_x;
    std::optional<Point_2> first_inter = first_intersection(ctx);

    if(auto y_at_txmin = ctx.bounded_compute_y_at_x(curve, txmin);
       y_at_txmin.has_value() && y_at_txmin != ctx->ymin() && y_at_txmin != ctx->ymax())
    {
      // The tight starting point of the curve is within the bbox and
      // it's not on the top or bottom edge.
      *ctx.out_it++ = txmin == ctx->xmin() ? ctx->approx_pt_on_boundary(Point_2(txmin, y_at_txmin.value()))
                                           : ctx->approx_pt(Point_2(txmin, y_at_txmin.value()));
      FT segment_end = first_inter.has_value() ? first_inter->x() : txmax;
      approximate_simple_curve_segment(ctx, txmin, segment_end, ctx->approx_error);
      last_x = segment_end;
    } else if(first_inter.has_value()) {
      last_x = first_inter->x();
    } else {
      // We assert that the curve is outbound.
      // If the min end is bounded, it's obvious.
      //
      // If the min end is unbounded, we know that the curve has no intersections with top, bottom or left edge(txmin ==
      // xmin when the min end is unbounded) and the min end of the curve is outbound (it approaches infinity in one or
      // both dimension).
      // Assume that the curve does has one point within the bbox. Note that it's a contiguous
      // x-monotone curve. So it must cross the top, bottom or left edge to reach the min end from right to left,
      // which is a contradiction.
      return polyline;
    }

    // iterate through the intersections and insert segments in-between.
    std::merge(top_inters.begin(), top_inters.end(), bottom_inters.begin(), bottom_inters.end(),
               boost::make_function_output_iterator([&last_x, &ctx](const Point_2& inter) {
                 approximate_simple_curve_segment(ctx, last_x, inter.x(), ctx->approx_error);
                 *ctx.out_it++ = ctx->approx_pt_on_boundary(inter);
                 last_x = inter.x();
               }),
               [&ctx](const Point_2& pt1, const Point_2& pt2) { return ctx->compare_xy_2(pt1, pt2) == CGAL::SMALLER; });

    if(auto y_at_txmax = ctx.bounded_compute_y_at_x(curve, txmax);
       y_at_txmax.has_value() && y_at_txmax != ctx->ymin() && y_at_txmax != ctx->ymax())
    {
      approximate_simple_curve_segment(ctx, last_x, txmax, ctx->approx_error);
      *ctx.out_it++ = txmax == ctx->xmax() ? ctx->approx_pt_on_boundary(Point_2(txmax, y_at_txmax.value()))
                                           : ctx->approx_pt(Point_2(txmax, y_at_txmax.value()));
    }

    return polyline;
  }

private:
  const Arr_bounded_render_context& m_ctx;
  const Arr_bounded_approximate_point_2& m_approx_pt;
  const Arr_bounded_compute_y_at_x m_bounded_compute_y_at_x;
  const X_monotone_curve_2 m_top;
  const X_monotone_curve_2 m_bottom;
};
} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_CURVE_2_H