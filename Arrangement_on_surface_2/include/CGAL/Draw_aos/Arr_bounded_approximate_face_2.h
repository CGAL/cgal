#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H

#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <algorithm>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_curve_2.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_point_2.h>
#include <CGAL/Draw_aos/Arr_bounded_face_triangulator.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {

namespace draw_aos {

/**
 * @brief A stateful geometry simplifier that simplifies horizontal and vertical segments
 *
 * @tparam GeomTraits
 * @tparam OutputIterator
 */
template <typename GeomTraits, typename OutputIterator>
class Colinear_simplifier
{
  using Geom_traits = GeomTraits;
  using Approx_point = typename Arr_approximate_traits<Geom_traits>::Approx_point;

public:
  using Insert_iterator = boost::function_output_iterator<std::function<void(Approx_point)>>;

public:
  Colinear_simplifier(OutputIterator out_it, Bbox_2 bbox)
      : m_out_it(out_it)
      , m_bbox(bbox) {}

  void dump() {
    if(m_start.has_value()) {
      *m_out_it++ = m_start.value();
      m_start.reset();
    }
    if(m_mid.has_value()) {
      *m_out_it++ = m_mid.value();
      m_mid.reset();
    }
  }

  Insert_iterator insert_iterator() {
    return boost::make_function_output_iterator(std::function([this](Approx_point p) {
      if(m_mid.has_value()) {
        if(p.y() == m_mid->y() && p.y() == m_start->y() || p.x() == m_mid->x() && p.x() == m_start->x())
          // Three points are collinear horizontally or vertically.
          m_mid = p;
        else {
          *m_out_it++ = m_start.value();
          m_start = m_mid;
          m_mid = p;
        }
        return;
      }

      if(m_start.has_value())
        m_mid = p;
      else
        m_start = p;
    }));
  }

  ~Colinear_simplifier() { dump(); }

private:
  OutputIterator m_out_it;
  std::optional<Approx_point> m_start, m_mid;
  const Bbox_2 m_bbox;
};

/**
 * @brief Bounded face approximation for arrangements.
 * @note Member functions are not thread-safe.
 */
template <typename Arrangement>
class Arr_bounded_approximate_face_2
{
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Polyline_geom = typename Approx_traits::Polyline_geom;
  using Ccb_halfedge_const_circulator = typename Arrangement::Ccb_halfedge_const_circulator;
  using Approx_point = typename Approx_traits::Approx_point;
  using Approx_nt = typename Approx_traits::Approx_nt;
  using Approx_kernel = typename Approx_traits::Approx_kernel;
  using Approx_line_2 = typename Approx_kernel::Line_2;
  using Triangulated_face = typename Approx_traits::Triangulated_face;
  using Feature_portal_map = typename Arr_portals<Arrangement>::Feature_portals_map;
  using Portal_exit = typename Arr_portals<Arrangement>::Portal_exit;
  using Portal_exit_vector = typename Arr_portals<Arrangement>::Portal_exit_vector;

  using Bounded_approximate_point_2 = Arr_bounded_approximate_point_2<Arrangement>;
  using Bounded_approximate_curve_2 = Arr_bounded_approximate_curve_2<Arrangement>;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

  using Triangulator = Arr_bounded_face_triangulator<Arrangement>;
  using Simplifier = Colinear_simplifier<Geom_traits, typename Triangulator::Insert_iterator>;

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
    using Output_iterator =
        typename Colinear_simplifier<Geom_traits, typename Triangulator::Insert_iterator>::Insert_iterator;

  public:
    using Insert_iterator = boost::function_output_iterator<std::function<void(Approx_point pt)>>;

    Context(const Bounded_render_context& ctx,
            Output_iterator out_it,
            const Bounded_approximate_point_2& bounded_approx_pt,
            const Bounded_approximate_curve_2& bounded_approx_curve)
        : Bounded_render_context(ctx)
        , m_base_out_it(out_it)
        , m_approx(ctx.m_traits.approximate_2_object())
        , m_bounded_approx_pt(bounded_approx_pt)
        , m_bounded_approx_curve(bounded_approx_curve) {
      m_out_it = boost::make_function_output_iterator(std::function([this](Approx_point pt) {
        if(!this->contains_x(pt.x())) return;
        pt = Approx_point(pt.x(), std::clamp(pt.y(), this->ymin(), this->ymax()));
        m_last_pt = pt;
        *m_base_out_it++ = pt;
      }));
    }
    // Let's not accidentally copy this object.
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

  private:
    Output_iterator m_base_out_it;

  public:
    const Bounded_approximate_point_2& m_bounded_approx_pt;
    const Bounded_approximate_curve_2& m_bounded_approx_curve;
    std::optional<Approx_point> m_last_pt;
    Insert_iterator m_out_it;
    const typename Geom_traits::Approximate_2 m_approx;
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
      Approx_point from_pt(ctx.m_last_pt.has_value() ? ctx.m_last_pt->x() : ctx.xmin(), ctx.ymin());
      Approx_point to_pt(ctx.xmax(), ctx.ymin());
      return Polyline_geom{from_pt, to_pt};
    }
    if(side == ARR_TOP_BOUNDARY) {
      Approx_point from_pt(ctx.xmin(), ctx.ymax());
      Approx_point to_pt(ctx.m_last_pt.has_value() ? ctx.m_last_pt->x() : ctx.xmax(), ctx.ymax());
      return Polyline_geom{from_pt, to_pt};
    }
    CGAL_assertion(false && "Unexpected side for a fictitious edge.");
    return Polyline_geom{};
  }

  static void approximate_vertex(Context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_at_open_boundary()) return;
    ctx.m_bounded_approx_pt(vh);
    // Handle portal from this vertex if there is one.
    auto portals_it = ctx.m_feature_portals.find(vh);
    if(portals_it == ctx.m_feature_portals.end()) return;
    const auto& portals = portals_it->second;
    if(portals.empty()) return;
    CGAL_assertion_msg(portals.size() == 1, "Vertex should have at most one portal.");
    traverse_portal(ctx, ctx.m_approx(vh->point()), portals.front());
  }

  template <typename DirectionTag>
  static void approximate_halfedge(Context& ctx, const Halfedge_const_handle& he, DirectionTag dir_tag) {
    constexpr bool Is_left_to_right = std::is_same_v<DirectionTag, Left_to_right_tag>;
    using Approx_point_it = std::conditional_t<Is_left_to_right, typename Polyline_geom::const_iterator,
                                               typename Polyline_geom::const_reverse_iterator>;
    using Portals_it = std::conditional_t<Is_left_to_right, typename Portal_exit_vector::const_iterator,
                                          typename Portal_exit_vector::const_reverse_iterator>;
    using Point_or_portal = std::variant<Approx_point, Vertex_const_handle>;

    const Polyline_geom& polyline =
        he->is_fictitious() ? approximate_fictitious_edge(ctx, he) : ctx.m_bounded_approx_curve(he);
    auto portal_map_it = ctx.m_feature_portals.find(he);
    const Portal_exit_vector& portals =
        portal_map_it != ctx.m_feature_portals.end() ? portal_map_it->second : Portal_exit_vector{};

    Approx_point_it points_iter, points_end;
    Portals_it portals_iter, portals_end;
    if constexpr(Is_left_to_right) {
      points_iter = polyline.begin();
      points_end = polyline.end();
      portals_iter = portals.begin();
      portals_end = portals.end();
    } else {
      points_iter = polyline.rbegin();
      points_end = polyline.rend();
      portals_iter = portals.rbegin();
      portals_end = portals.rend();
    }
    auto less_than = [](Approx_nt x1, Approx_nt x2) { return Is_left_to_right ? x1 < x2 : x1 > x2; };

    // A variant of two pointers algorithm, but std::merge() can't fit in our purpose.
    std::optional<Approx_point> prev_pt;
    for(; points_iter != points_end; ++points_iter) {
      Approx_point curr_pt = *points_iter;

      for(; portals_iter != portals_end; ++portals_iter) {
        const Portal_exit& exit = *portals_iter;
        Approx_point exit_pt = ctx.m_approx(exit->point());
        if(!less_than(exit_pt.x(), curr_pt.x()) || !prev_pt.has_value()) break;
        if(less_than(exit_pt.x(), prev_pt->x())) {
          traverse_portal(ctx, exit_pt, exit);
          continue;
        }
        // Compute the portal entry point.
        std::optional<std::variant<Approx_point, Approx_line_2>> res = CGAL::intersection(
            Approx_line_2(*prev_pt, curr_pt), Approx_line_2(exit_pt, Approx_point(exit_pt.x(), exit_pt.y() + 1)));
        CGAL_assertion(res.has_value() && std::holds_alternative<Approx_point>(*res));
        traverse_portal(ctx, std::get<Approx_point>(*res), exit);
      }

      *ctx.m_out_it++ = curr_pt;
      prev_pt = curr_pt;
    }

    for(; portals_iter != portals_end; ++portals_iter) {
      const Portal_exit& exit = *portals_iter;
      traverse_portal(ctx, ctx.m_approx(exit->point()), exit);
    }
  }

  static void approximate_inner_ccb(Context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_isolated()) {
      approximate_vertex(ctx, vh);
      return;
    }
    approximate_ccb(ctx, get_inner_ccb_first_halfedge(vh));
  }

  static void traverse_portal(Context& ctx, Approx_point entry, const Portal_exit& exit) {
    *ctx.m_out_it++ = entry;
    approximate_inner_ccb(ctx, exit);
    // We come back here after inner ccb is approximated.
    *ctx.m_out_it++ = entry;
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
  Arr_bounded_approximate_face_2(const Bounded_render_context& ctx,
                                 const Bounded_approximate_point_2& bounded_approx_pt,
                                 const Bounded_approximate_curve_2& bounded_approx_curve)
      : m_ctx(ctx)
      , m_bounded_approx_pt(bounded_approx_pt)
      , m_bounded_approx_curve(bounded_approx_curve) {}

  const Triangulated_face& operator()(const Face_const_handle& fh) const {
    auto [triangulated_face, inserted] = m_ctx.m_cache.try_emplace(fh);
    if(!inserted) return triangulated_face;
    if(m_ctx.is_cancelled()) return triangulated_face;

    CGAL_precondition_msg(!fh->is_fictitious(), "Cannot approximate a fictitious face.");

    if(!fh->has_outer_ccb()) {
      // The face is the unbounded face of bounded arrangements, we skip approximating any non degenerate features.
      for(auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
        auto circ = *inner_ccb;
        do {
          if(circ->face() != circ->twin()->face()) continue;
          m_bounded_approx_curve(circ);
        } while(++circ != *inner_ccb);
      }
      for(auto vh = fh->isolated_vertices_begin(); vh != fh->isolated_vertices_end(); ++vh) m_bounded_approx_pt(vh);
      return triangulated_face;
    }

    auto triangulator = Triangulator(m_ctx);
    auto simplifier = Simplifier(triangulator.insert_iterator(), m_ctx.bbox());
    auto ctx = Context(m_ctx, simplifier.insert_iterator(), m_bounded_approx_pt, m_bounded_approx_curve);
    approximate_ccb(ctx, fh->outer_ccb());
    simplifier.dump();
    return triangulated_face = std::move(triangulator);
  }

private:
  const Bounded_render_context& m_ctx;
  const Bounded_approximate_point_2& m_bounded_approx_pt;
  const Bounded_approximate_curve_2& m_bounded_approx_curve;
};

} // namespace draw_aos
} // namespace CGAL

#endif