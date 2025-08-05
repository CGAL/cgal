#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_H

#include <cstddef>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <algorithm>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_halfedge.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_vertex.h>
#include <CGAL/Draw_aos/Arr_projector.h>
#include <CGAL/Draw_aos/Arr_bounded_face_triangulator.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>
#include "CGAL/basic.h"
#include "CGAL/enum.h"
#include "CGAL/number_utils.h"

namespace CGAL {

namespace draw_aos {

/**
 * @brief Bounded face approximation for arrangements.
 * @note Member functions are not thread-safe.
 */
template <typename Arrangement>
class Arr_bounded_approximate_face
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Point_geom = typename Approx_traits::Point_geom;
  using Polyline_geom = typename Approx_traits::Polyline_geom;
  using Point_2 = typename Geom_traits::Point_2;
  using Ccb_halfedge_const_circulator = typename Arrangement::Ccb_halfedge_const_circulator;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Approx_kernel = typename Approx_traits::Approx_kernel;
  using Approx_line_2 = typename Approx_kernel::Line_2;
  using Triangulated_face = typename Approx_traits::Triangulated_face;
  using Portal_exit = typename Arr_portals<Arrangement>::Portal_exit;
  using Portal_exit_vector = typename Arr_portals<Arrangement>::Portal_exit_vector;

  using Bounded_approximate_vertex = Arr_bounded_approximate_vertex<Arrangement>;
  using Bounded_approximate_halfedge = Arr_bounded_approximate_halfedge<Arrangement>;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

  using Triangulator = Arr_bounded_face_triangulator<Arrangement>;

  struct Left_to_right_tag
  {};
  struct Right_to_left_tag
  {};

  struct Inner_ccb_tag
  {};
  struct Outer_ccb_tag
  {};

private:
  class Context : public Bounded_render_context
  {
  private:
    using Output_iterator = typename Triangulator::Insert_iterator;

  public:
    using Insert_iterator = boost::function_output_iterator<std::function<void(Point_geom pt)>>;

    Context(const Bounded_render_context& ctx,
            Output_iterator out_it,
            const Bounded_approximate_vertex& bounded_approx_vertex,
            const Bounded_approximate_halfedge& bounded_approx_halfedge)
        : Bounded_render_context(ctx)
        , m_base_out_it(out_it)
        , m_approx(ctx.m_traits.approximate_2_object())
        , m_proj(ctx.m_traits)
        , m_bounded_approx_vertex(bounded_approx_vertex)
        , m_bounded_approx_halfedge(bounded_approx_halfedge) {
      m_out_it = boost::make_function_output_iterator(std::function([this](Point_geom pt) {
        if(!this->contains_x(pt.x())) return;
        pt = Point_geom(pt.x(), std::clamp(pt.y(), this->ymin(), this->ymax()));
        m_last_pt = pt;
        *m_base_out_it++ = pt;
      }));
    }
    // Let's not accidentally copy this object.
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

    Point_geom approx(const Point_2& pt) const { return m_proj.project(m_approx(pt)); }

  private:
    Output_iterator m_base_out_it;
    const typename Geom_traits::Approximate_2 m_approx;
    const Arr_projector<Geom_traits> m_proj;

  public:
    const Bounded_approximate_vertex& m_bounded_approx_vertex;
    const Bounded_approximate_halfedge& m_bounded_approx_halfedge;
    std::optional<Point_geom> m_last_pt;
    Insert_iterator m_out_it;
  };

private:
  /**
   * @brief Get the first halfedge on the inner ccb that should be traversed.
   *
   * @precondition: The vertex must not be isolated.
   * @param vh
   */
  static Halfedge_const_handle get_inner_ccb_first_halfedge(const Vertex_const_handle& vh) {
    auto circ = vh->incident_halfedges();
    if(vh->degree() == 1) return circ->twin();

    auto curr = circ;
    auto next = curr;
    ++next;
    // Traverse the halfedges incident to the vertex in clockwise direction.
    // Note that circulator around a vertex points to a halfedge whose target is the vertex.
    do {
      // If "curr" crosses 12 o'clock and reaches "next", we found the wanted halfedge.
      if(curr->direction() == ARR_LEFT_TO_RIGHT && next->direction() == ARR_RIGHT_TO_LEFT) break;
      ++curr;
      ++next;
    } while(curr != circ);
    // The opposite halfedge is the one we want.
    return next->twin();
  }

  static Arr_parameter_space side_of_fictitious_edge(const Halfedge_const_handle& he) {
    const auto& source = he->source();
    const auto& target = he->target();
    Arr_parameter_space src_x_space = source->parameter_space_in_x();
    Arr_parameter_space src_y_space = source->parameter_space_in_y();
    Arr_parameter_space tgt_x_space = target->parameter_space_in_x();
    Arr_parameter_space tgt_y_space = target->parameter_space_in_y();
    if(src_x_space == tgt_x_space && src_x_space != ARR_INTERIOR) return src_x_space;
    if(src_y_space == tgt_y_space && src_y_space != ARR_INTERIOR) return src_y_space;
    CGAL_assertion(false && "Unexpected parameter space for fictitious edge vertices.");
    return ARR_INTERIOR;
  }

  // Generate dummy curve(directed left to right) for the fictitious edge.
  static Polyline_geom approximate_fictitious_edge(const Context& ctx, const Halfedge_const_handle& he) {
    Arr_parameter_space side = side_of_fictitious_edge(he);
    // There's no need to handle fictitious edges on left or right boundaries.
    if(side == ARR_LEFT_BOUNDARY || side == ARR_RIGHT_BOUNDARY) return Polyline_geom{};
    if(side == ARR_BOTTOM_BOUNDARY) {
      Point_geom from_pt(ctx.m_last_pt.has_value() ? ctx.m_last_pt->x() : ctx.xmin(), ctx.ymin());
      Point_geom to_pt(ctx.xmax(), ctx.ymin());
      return Polyline_geom{from_pt, to_pt};
    }
    if(side == ARR_TOP_BOUNDARY) {
      Point_geom from_pt(ctx.xmin(), ctx.ymax());
      Point_geom to_pt(ctx.m_last_pt.has_value() ? ctx.m_last_pt->x() : ctx.xmax(), ctx.ymax());
      return Polyline_geom{from_pt, to_pt};
    }
    CGAL_assertion(false && "Unexpected side for a fictitious edge.");
    return Polyline_geom{};
  }

  static void approximate_vertex(Context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_at_open_boundary()) return;
    ctx.m_bounded_approx_vertex(vh);
    // Handle portal from this vertex if there is one.
    auto portals_it = ctx.m_portals_map.find(vh);
    if(portals_it == ctx.m_portals_map.end()) return;
    const auto& portals = portals_it->second;
    if(portals.empty()) return;
    CGAL_assertion_msg(portals.size() == 1, "Vertex should have at most one portal.");
    traverse_portal(ctx, ctx.approx(vh->point()), portals.front());
  }

  template <typename DirectionTag>
  static void approximate_halfedge(Context& ctx, const Halfedge_const_handle& he, DirectionTag dir_tag) {
    const Polyline_geom& polyline =
        he->is_fictitious() ? approximate_fictitious_edge(ctx, he) : ctx.m_bounded_approx_halfedge(he);
    auto portal_map_it = ctx.m_portals_map.find(he);
    const Portal_exit_vector& portals =
        portal_map_it != ctx.m_portals_map.end() ? portal_map_it->second : Portal_exit_vector{};

    // A variant of two pointers algorithm. std::merge() can't fit in our purpose so we implement it ourselves.
    auto compare_x = [](const Point_geom& a, const Point_geom& b) -> Comparison_result {
      return std::is_same_v<DirectionTag, Left_to_right_tag> ? sign(a.x() - b.x()) : sign(b.x() - a.x());
    };
    std::optional<Point_geom> prev_pt;
    auto portals_iter = portals.begin(), portals_end = portals.end();
    for(const auto& curr_pt : polyline) {
      for(; portals_iter != portals_end; ++portals_iter) {
        const Portal_exit& exit = *portals_iter;
        Point_geom exit_pt = ctx.approx(exit->point());
        if(compare_x(exit_pt, curr_pt) >= 0 || !prev_pt.has_value()) break;
        if(compare_x(exit_pt, *prev_pt) < 0) {
          traverse_portal(ctx, exit_pt, exit);
          continue;
        }
        // Compute the portal entry point.
        std::optional<std::variant<Point_geom, Approx_line_2>> res = CGAL::intersection(
            Approx_line_2(*prev_pt, curr_pt), Approx_line_2(exit_pt, Point_geom(exit_pt.x(), exit_pt.y() + 1)));
        CGAL_assertion(res.has_value() && std::holds_alternative<Point_geom>(*res));
        traverse_portal(ctx, std::get<Point_geom>(*res), exit);
      }

      *ctx.m_out_it++ = curr_pt;
      prev_pt = curr_pt;
    }

    for(; portals_iter != portals_end; ++portals_iter) {
      const Portal_exit& exit = *portals_iter;
      traverse_portal(ctx, ctx.approx(exit->point()), exit);
    }
  }

  static void approximate_inner_ccb(Context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_isolated()) {
      approximate_vertex(ctx, vh);
      return;
    }
    approximate_ccb(ctx, get_inner_ccb_first_halfedge(vh));
  }

  static void traverse_portal(Context& ctx, Point_geom entry_pt, const Vertex_const_handle& exit_vh) {
    Point_geom exit_pt = ctx.approx(exit_vh->point());
    Approx_nt portal_x = exit_pt.x();
    // Try enforce a vertical segment.
    entry_pt = Point_geom(portal_x, entry_pt.y());
    *ctx.m_out_it++ = entry_pt;

    if constexpr(is_or_derived_from_curved_surf_traits<Geom_traits>)
      // Approximate the vertical segment, only meaningful when dealing with arranagements on curved surfaces.
      for(Approx_nt y = entry_pt.y() - ctx.m_approx_error; y > exit_pt.y(); y -= ctx.m_approx_error)
        *ctx.m_out_it++ = Point_geom(portal_x, y);

    approximate_inner_ccb(ctx, exit_vh);
    // Come back here after inner ccb is approximated.
    *ctx.m_out_it++ = entry_pt;
  }

  static void approximate_ccb(Context& ctx, Ccb_halfedge_const_circulator start_circ) {
    // Try to start on a concrete halfedge.
    // For any unbounded face, there can't be more than 4 adjacent fictitious edges.
    for(int i = 0; i < 4 && start_circ->is_fictitious(); ++i) ++start_circ;

    auto circ = start_circ;
    do {
      circ->direction() == ARR_LEFT_TO_RIGHT ? approximate_halfedge(ctx, circ, Left_to_right_tag{})
                                             : approximate_halfedge(ctx, circ, Right_to_left_tag{});
      approximate_vertex(ctx, circ->target());
    } while(++circ != start_circ);
  }

public:
  Arr_bounded_approximate_face(const Bounded_render_context& ctx)
      : m_ctx(ctx)
      , m_bounded_approx_vertex(ctx)
      , m_bounded_approx_halfedge(ctx) {}

  const Triangulated_face& operator()(const Face_const_handle& fh) const {
    CGAL_precondition_msg(!fh->is_fictitious(), "Cannot approximate a fictitious face.");

    auto [triangulated_face, inserted] = m_ctx.m_cache.try_emplace(fh);
    if(!inserted) return triangulated_face;
    if(m_ctx.is_cancelled()) return triangulated_face;

    if(!fh->has_outer_ccb()) {
      // The face is the unbounded face of bounded arrangements, we skip approximating any non degenerate features.
      for(auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
        auto circ = *inner_ccb;
        do {
          if(circ->face() != circ->twin()->face()) continue;
          m_bounded_approx_halfedge(circ);
        } while(++circ != *inner_ccb);
      }
      for(auto vh = fh->isolated_vertices_begin(); vh != fh->isolated_vertices_end(); ++vh) m_bounded_approx_vertex(vh);
      return triangulated_face;
    }

    auto triangulator = Triangulator(m_ctx);
    auto ctx = Context(m_ctx, triangulator.insert_iterator(), m_bounded_approx_vertex, m_bounded_approx_halfedge);
    approximate_ccb(ctx, fh->outer_ccb());
    return triangulated_face = std::move(triangulator);
  }

private:
  const Bounded_render_context& m_ctx;
  const Bounded_approximate_vertex m_bounded_approx_vertex;
  const Bounded_approximate_halfedge m_bounded_approx_halfedge;
};

} // namespace draw_aos
} // namespace CGAL

#endif