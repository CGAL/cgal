#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_CURVE_2_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_CURVE_2_H

#include <array>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <optional>

#include <CGAL/Arr_enums.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>

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
class Arr_bounded_approximate_curve_2
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;

  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approx_kernel = typename Approx_traits::Approx_kernel;
  using Approx_line_2 = typename Approx_kernel::Line_2;
  using Polyline_geom = typename Approx_traits::Polyline_geom;

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
        , m_approx(ctx.m_traits.approximate_2_object())
        , m_base_out_it(std::back_inserter(polyline))
        , m_out_it(boost::make_function_output_iterator(std::function([this](Approx_point pt) {
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
    const X_monotone_curve_2& m_curve;
    const typename Geom_traits::Approximate_2 m_approx;
    const Boundary_lines& m_boundary_lines;
    std::optional<Approx_point> m_last_pt;
    boost::function_output_iterator<std::function<void(Approx_point)>> m_out_it;
  };

  static Approx_point trace_boundary_inter(const Context& ctx, Approx_point pt, Side_of_boundary side) {
    Approx_point inter = std::get<Approx_point>(
        *CGAL::intersection(Approx_line_2(*ctx.m_last_pt, pt), ctx.m_boundary_lines[static_cast<std::size_t>(side)]));
    // Prevent floating point errors.
    switch(side) {
    case Side_of_boundary::Left:
      return Approx_point(ctx.xmin(), inter.y());
    case Side_of_boundary::Right:
      return Approx_point(ctx.xmax(), inter.y());
    case Side_of_boundary::Top:
      return Approx_point(inter.x(), ctx.ymax());
    case Side_of_boundary::Bottom:
      return Approx_point(inter.x(), ctx.ymin());
    default:
      CGAL_assertion(false && "Unexpected side of boundary.");
      return Approx_point();
    }
  }

  static void update(Context& ctx, Approx_point pt) {
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
  Arr_bounded_approximate_curve_2(const Bounded_render_context& ctx)
      : m_ctx(ctx)
      , m_boundary_lines({
            Approx_line_2(Approx_point(ctx.xmin(), ctx.ymax()), Approx_point(ctx.xmax(), ctx.ymax())), // Top = 0
            Approx_line_2(Approx_point(ctx.xmin(), ctx.ymin()), Approx_point(ctx.xmin(), ctx.ymax())), // Left = 1
            Approx_line_2(Approx_point(ctx.xmin(), ctx.ymin()), Approx_point(ctx.xmax(), ctx.ymin())), // Bottom = 2
            Approx_line_2(Approx_point(ctx.xmax(), ctx.ymin()), Approx_point(ctx.xmax(), ctx.ymax())), // Right = 3
        }) {}

  const Polyline_geom& operator()(const Halfedge_const_handle& he) const {
    CGAL_assertion(!he->is_fictitious());

    auto [polyline, inserted] = m_ctx.m_cache.try_emplace(he);
    if(!inserted) return polyline;
    if(m_ctx.is_cancelled()) return polyline;

    const X_monotone_curve_2& curve = he->curve();
    Context ctx(m_ctx, curve, polyline, m_boundary_lines);
    m_ctx.m_traits.approximate_2_object()(
        curve, m_ctx.m_approx_error, boost::make_function_output_iterator([&ctx](Approx_point pt) { update(ctx, pt); }),
        true);
    return polyline;
  }

private:
  const Bounded_render_context& m_ctx;
  const Boundary_lines m_boundary_lines;
};

} // namespace draw_aos
} // namespace CGAL
#endif