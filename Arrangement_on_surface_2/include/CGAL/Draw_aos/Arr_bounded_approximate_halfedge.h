#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_HALFEDGE_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_HALFEDGE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <optional>
#include <type_traits>

#include <boost/range/any_range.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>
#include "CGAL/Draw_aos/Arr_projector.h"

namespace CGAL {
namespace draw_aos {

/**
 * @brief Functor to approximate an x-monotone curve within an bounding box.
 *
 * The Approximation is done from xmin to xmax with a given step. For parts outbound the y limits and precedes or
 * succeeds a part within, the approximation may be skipped but there will be at least one point outside the bbox
 * for indication.
 */
template <typename Arrangement>
class Arr_bounded_approximate_halfedge
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;

  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Approx_point = typename Approx_traits::Approx_point;
  using Point_geom = typename Approx_traits::Point_geom;
  using Polyline_geom = typename Approx_traits::Polyline_geom;
  using Approx_kernel = typename Approx_traits::Approx_kernel;
  using Approx_line_2 = typename Approx_kernel::Line_2;

  using X_monotone_curve_2 = typename Geom_traits::X_monotone_curve_2;
  using Point_2 = typename Geom_traits::Point_2;

  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;
  using Boundary_lines = std::array<Approx_line_2, 4>;

private:
  struct Context : public Bounded_render_context
  {
    Context(const Bounded_render_context& ctx,
            const X_monotone_curve_2& curve,
            Polyline_geom& polyline,
            const Boundary_lines& boundary_lines)
        : Bounded_render_context(ctx)
        , m_curve(curve)
        , m_boundary_lines(boundary_lines)
        , m_proj(ctx.m_traits)
        , m_base_out_it(std::back_inserter(polyline))
        , m_out_it(boost::make_function_output_iterator(std::function([this](Point_geom pt) {
          if(pt.x() < this->xmin()) {
            // We need the last point if not yet x-inbound.
            m_last_pt = pt;
            return;
          } else if(pt.x() > this->xmax())
            return;

          *m_base_out_it++ = pt;
          m_last_pt = pt;
        }))) {}

  private:
    std::back_insert_iterator<Polyline_geom> m_base_out_it;

  public:
    Arr_projector<Geom_traits> m_proj;
    const X_monotone_curve_2& m_curve;
    const Boundary_lines& m_boundary_lines;
    std::optional<Point_geom> m_last_pt;
    boost::function_output_iterator<std::function<void(Point_geom)>> m_out_it;
  };

  static Point_geom trace_boundary_inter(const Context& ctx, Point_geom pt, Side_of_boundary side) {
    Point_geom inter = std::get<Point_geom>(
        *CGAL::intersection(Approx_line_2(*ctx.m_last_pt, pt), ctx.m_boundary_lines[static_cast<std::size_t>(side)]));
    // Prevent floating point errors.
    switch(side) {
    case Side_of_boundary::Left:
      return Point_geom(ctx.xmin(), inter.y());
    case Side_of_boundary::Right:
      return Point_geom(ctx.xmax(), inter.y());
    case Side_of_boundary::Top:
      return Point_geom(inter.x(), ctx.ymax());
    case Side_of_boundary::Bottom:
      return Point_geom(inter.x(), ctx.ymin());
    default:
      CGAL_assertion(false && "Unexpected side of boundary.");
      return Point_geom();
    }
  }

  static void update(Context& ctx, Point_geom pt) {
    if(!ctx.m_last_pt.has_value()) {
      *ctx.m_out_it++ = pt;
      return;
    }

    if(ctx.m_last_pt->x() < ctx.xmin() && pt.x() >= ctx.xmin()) {
      *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Left);
    }

    if(ctx.m_last_pt->y() < ctx.ymin()) {
      if(pt.y() > ctx.ymin()) {
        *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Bottom);
      }
      if(pt.y() > ctx.ymax()) {
        *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Top);
      }
    } else if(ctx.m_last_pt->y() > ctx.ymax()) {
      if(pt.y() < ctx.ymax()) {
        *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Top);
      }
      if(pt.y() < ctx.ymin()) {
        *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Bottom);
      }
    } else {
      if(pt.y() < ctx.ymin()) {
        *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Bottom);
      } else if(pt.y() > ctx.ymax()) {
        *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Top);
      }
    }

    if(ctx.m_last_pt->x() <= ctx.xmax() && pt.x() > ctx.xmax()) {
      *ctx.m_out_it++ = trace_boundary_inter(ctx, pt, Side_of_boundary::Right);
    }

    *ctx.m_out_it++ = pt;
  }

public:
  Arr_bounded_approximate_halfedge(const Bounded_render_context& ctx)
      : m_ctx(ctx)
      , m_boundary_lines({
            Approx_line_2(Point_geom(ctx.xmin(), ctx.ymax()), Point_geom(ctx.xmax(), ctx.ymax())), // Top = 0
            Approx_line_2(Point_geom(ctx.xmin(), ctx.ymin()), Point_geom(ctx.xmin(), ctx.ymax())), // Left = 1
            Approx_line_2(Point_geom(ctx.xmin(), ctx.ymin()), Point_geom(ctx.xmax(), ctx.ymin())), // Bottom = 2
            Approx_line_2(Point_geom(ctx.xmax(), ctx.ymin()), Point_geom(ctx.xmax(), ctx.ymax())), // Right = 3
        }) {}

  const Polyline_geom& operator()(const Halfedge_const_handle& he) const {
    CGAL_assertion(!he->is_fictitious());

    auto [polyline, inserted] = m_ctx.m_cache.try_emplace(he);
    if(!inserted) return polyline;
    if(m_ctx.is_cancelled()) return polyline;

    const X_monotone_curve_2& curve = he->curve();
    Context ctx(m_ctx, curve, polyline, m_boundary_lines);
    m_ctx.m_traits.approximate_2_object()(
        curve, m_ctx.m_approx_error,
        boost::make_function_output_iterator([&ctx](Approx_point pt) { update(ctx, ctx.m_proj.project(pt)); }),
        true); // ltr ordering

    // also approximate the twin halfedge
    auto [twin_poly, twin_inserted] = m_ctx.m_cache.try_emplace(he->twin());
    twin_poly = polyline;
    if(twin_inserted) adapt_polyline(ctx, twin_poly, he->twin());

    adapt_polyline(ctx, polyline, he);
    return polyline;
  }

  /**
   * @brief Functor to adapt approximated curve points based on the actual halfedge, giving correct ordering,
   * continuity, etc.
   */
  static void adapt_polyline(Context& ctx, Polyline_geom& polyline, const Halfedge_const_handle& he) {
    adapt_polyline_impl<Geom_traits>(ctx, polyline, he);
  }

  template <typename T, std::enable_if_t<!is_or_derived_from_agas_v<T>, int> = 0>
  static void adapt_polyline_impl(Context& ctx, Polyline_geom& polyline, const Halfedge_const_handle& he) {
    if(he->direction() == CGAL::ARR_LEFT_TO_RIGHT) return;
    std::reverse(polyline.begin(), polyline.end());
  }

  template <typename T, std::enable_if_t<is_or_derived_from_agas_v<T>, int> = 0>
  static void adapt_polyline_impl(Context& ctx, Polyline_geom& polyline, const Halfedge_const_handle& he) {
    if(polyline.size() < 2) return;

    auto azimuth_of_vertical_curve = [&ctx](const X_monotone_curve_2& curve) -> Approx_nt {
      using Point_2 = typename Geom_traits::Point_2;
      using Direction_3 = typename Geom_traits::Direction_3;

      const auto& traits = ctx.m_traits;
      const Direction_3& normal_dir = curve.normal();
      if(normal_dir.dx() == 0 && normal_dir.dy() == 1) return 0; // overlaps with the identification curve
      Point_2 normal_pt(normal_dir, Point_2::NO_BOUNDARY_LOC);
      Approx_point normal_approx_pt = traits.approximate_2_object()(normal_pt);
      Point_geom normal_proj_pt = ctx.m_proj.project(normal_approx_pt);
      return std::fmod(normal_proj_pt.x() + 0.5 * CGAL_PI, CGAL_PI * 2.0);
    };

    const X_monotone_curve_2& curve = he->curve();

    if(curve.is_vertical()) {
      Approx_nt azimuth = azimuth_of_vertical_curve(curve);
      if(azimuth == 0 && he->direction() == ARR_LEFT_TO_RIGHT) azimuth = 2 * CGAL_PI;
      std::transform(polyline.begin(), polyline.end(), polyline.begin(),
                     [azimuth](Point_geom pt) { return Point_geom(azimuth, pt.y()); });
    } else if(polyline.back().x() == 0) {
      // For strictly x-monotone arcs, if the target point sits on the boundary, the x should be set to 2 * CGAL_PI
      polyline.back() = Point_geom(2 * CGAL_PI, polyline.back().y());
    }

    if(he->direction() == ARR_RIGHT_TO_LEFT) std::reverse(polyline.begin(), polyline.end());
  }

private:
  const Bounded_render_context& m_ctx;
  const Boundary_lines m_boundary_lines;
};

} // namespace draw_aos
} // namespace CGAL
#endif