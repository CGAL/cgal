#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H

#include <cstddef>
#include <functional>
#include <optional>
#include <type_traits>
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
 * @brief Patches corners between two boundary points of the bbox
 * counter-clockwisely.
 */
template <typename GeomTraits>
class Patch_boundary
{
  using Approx_point = typename Arr_approximation_geometry_traits<GeomTraits>::Approx_point;

private:
  Side_of_boundary side_of_boundary(Approx_point pt) const {
    if(m_bbox.xmin() <= pt.x() && pt.x() < m_bbox.xmax() && pt.y() == m_bbox.ymax()) {
      return Side_of_boundary::Top;
    } else if(pt.x() == m_bbox.xmin() && m_bbox.ymin() <= pt.y() && pt.y() < m_bbox.ymax()) {
      return Side_of_boundary::Left;
    } else if(m_bbox.xmin() < pt.x() && pt.x() <= m_bbox.xmax() && pt.y() == m_bbox.ymin()) {
      return Side_of_boundary::Bottom;
    } else if(pt.x() == m_bbox.xmax() && m_bbox.ymin() < pt.y() && pt.y() <= m_bbox.ymax()) {
      return Side_of_boundary::Right;
    } else {
      return Side_of_boundary::None;
    }
  }

  Approx_point corner_of_side(Side_of_boundary side) const {
    switch(side) {
    case Side_of_boundary::Top:
      // return the top-left corner
      return Approx_point(m_bbox.xmin(), m_bbox.ymax());
    case Side_of_boundary::Left:
      // return the bottom-left corner
      return Approx_point(m_bbox.xmin(), m_bbox.ymin());
    case Side_of_boundary::Bottom:
      // return the bottom-right corner
      return Approx_point(m_bbox.xmax(), m_bbox.ymin());
    case Side_of_boundary::Right:
      // return the top-right corner
      return Approx_point(m_bbox.xmax(), m_bbox.ymax());
    default:
      CGAL_assertion(false && "Invalid side of rectangle");
      return Approx_point();
    }
  }

  Side_of_boundary prev_side(Side_of_boundary side) const {
    CGAL_assertion(side != Side_of_boundary::None);
    return static_cast<Side_of_boundary>((static_cast<int>(side) + 3) % 4);
  }

  Side_of_boundary next_side(Side_of_boundary side) const {
    CGAL_assertion(side != Side_of_boundary::None);
    return static_cast<Side_of_boundary>((static_cast<int>(side) + 1) % 4);
  }

  // Computes the distance between two points with the precondition that they are on the same side of the boundary.
  double distance_on_same_side(Approx_point pt1, Approx_point pt2) const {
    return std::abs(pt1.x() - pt2.x()) + std::abs(pt1.y() - pt2.y());
  }

public:
  Patch_boundary(Bbox_2 bbox)
      : m_bbox(bbox) {}

  /**
   * @brief Patch the boundary between two points on the boundary of the bbox counter-clockwisely.
   */
  template <typename OutputIterator>
  void operator()(Approx_point from, Approx_point to, OutputIterator out_it) const {
    auto from_side = side_of_boundary(from);
    auto to_side = side_of_boundary(to);

    if(from_side == Side_of_boundary::None || to_side == Side_of_boundary::None) {
      return;
    }

    bool from_is_before_to = from_side == to_side && distance_on_same_side(from, corner_of_side(from_side)) <=
                                                         distance_on_same_side(to, corner_of_side(to_side));
    // Special case: if both on the same side and the current point is before
    // the last point in counter-clockwise boundary of the bbox, We have to go
    // around the bbox. Note that this includes when two points overlaps
    int num_corners_to_patch =
        from_is_before_to ? 4 : (static_cast<int>(to_side) - static_cast<int>(from_side) + 4) % 4;

    for(int i = 0; i < num_corners_to_patch; ++i) {
      Side_of_boundary side = static_cast<Side_of_boundary>((static_cast<int>(from_side) + i) % 4);
      Approx_point corner = corner_of_side(side);
      if(corner == from) {
        continue;
      }
      *out_it++ = corner;
    }
    return;
  }

  /**
   * @brief Patch all four sides.
   */
  template <typename OutputIterator>
  void operator()(OutputIterator out_it) const {
    for(auto side : {Side_of_boundary::Top, Side_of_boundary::Left, Side_of_boundary::Bottom, Side_of_boundary::Right})
    {
      *out_it++ = corner_of_side(side);
    }
  }

private:
  const Bbox_2 m_bbox;
};

template <typename GeomTraits, typename OutputIterator>
class Colinear_simplifier
{
  using Geom_traits = GeomTraits;
  using Approx_point = typename Arr_approximation_geometry_traits<Geom_traits>::Approx_point;

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
        if(p.y() == m_mid->y() && p.y() == m_start->y() || p.x() == m_mid->x() && p.x() == m_start->x()) {
          // Three points are collinear horizontally or vertically.
          m_mid = p;
        } else {
          *m_out_it++ = m_start.value();
          m_start = m_mid;
          m_mid = p;
        }
        return;
      }

      if(m_start.has_value()) {
        m_mid = p;
      } else {
        m_start = p;
      }
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
  using Approx_traits = Arr_approximation_geometry_traits<Geom_traits>;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Polyline_geom = typename Approx_traits::Polyline_geom;
  using Ccb_halfedge_const_circulator = typename Arrangement::Ccb_halfedge_const_circulator;
  using Approx_point = typename Approx_traits::Approx_point;
  using Triangulated_face = typename Approx_traits::Triangulated_face;
  using Feature_portal_map = typename Arr_portals<Arrangement>::Feature_portals_map;
  using Portal_vector = typename Arr_portals<Arrangement>::Portal_vector;
  using Portal = typename Arr_portals<Arrangement>::Portal;
  using Point_or_portal = std::variant<Approx_point, Portal>;
  using FT = typename Traits_adaptor<Geom_traits>::FT;

  using Bounded_approximate_point_2 = Arr_bounded_approximate_point_2<Arrangement>;
  using Bounded_approximate_curve_2 = Arr_bounded_approximate_curve_2<Arrangement>;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

  using Triangulator = Arr_bounded_face_triangulator<Arrangement>;
  using Patch = Patch_boundary<Geom_traits>;
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
  class Execution_context : public Arr_context_delegator<Bounded_render_context>
  {
  private:
    using Output_iterator =
        typename Colinear_simplifier<Geom_traits, typename Triangulator::Insert_iterator>::Insert_iterator;

  public:
    using Insert_iterator = boost::function_output_iterator<std::function<void(Approx_point pt)>>;

    Execution_context(const Bounded_render_context& ctx,
                      const Patch& patch,
                      Output_iterator out_it,
                      const Bounded_approximate_point_2& bounded_approx_pt,
                      const Bounded_approximate_curve_2& bounded_approx_curve)
        : Arr_context_delegator<Bounded_render_context>(ctx)
        , base_out_it(out_it)
        , patch(patch)
        , bounded_approx_pt(bounded_approx_pt)
        , bounded_approx_curve(bounded_approx_curve) {
      this->out_it = boost::make_function_output_iterator(std::function([this](Approx_point pt) {
        if(!(*this)->contains_x(pt.x())) {
          return;
        }

        pt = Approx_point(pt.x(), std::clamp(pt.y(), (*this)->ymin(), (*this)->ymax()));
        if(pt == last_pt) {
          passed_fictitious_edge = false;
          return;
        }

        if(passed_fictitious_edge && last_pt.has_value()) {
          this->patch(last_pt.value(), pt, base_out_it);
        }
        passed_fictitious_edge = false;

        *base_out_it++ = pt;

        if(!first_pt.has_value()) {
          first_pt = pt;
        }
        last_pt = pt;
      }));
    }
    // Let's not accidentally copy this object.
    Execution_context(const Execution_context&) = delete;
    Execution_context& operator=(const Execution_context&) = delete;

  public:
    const Bounded_approximate_point_2& bounded_approx_pt;
    const Bounded_approximate_curve_2& bounded_approx_curve;
    Insert_iterator out_it;
    std::optional<Approx_point> last_pt, first_pt;
    bool passed_fictitious_edge{false};
    Output_iterator base_out_it;
    const Patch& patch;
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
    if(vh->degree() == 1) {
      return circ->twin();
    }
    auto curr = circ;
    auto next = curr;
    ++next;
    // Traverse the halfedges incident to the vertex in clockwise direction.
    // Note that circulator around a vertex points to a halfedge whose target is the vertex.
    do {
      if(curr->direction() == ARR_LEFT_TO_RIGHT && next->direction() == ARR_RIGHT_TO_LEFT) {
        // This indicates that "curr" crosses 12 o'clock and reaches "next".
        break;
      }
      ++curr;
      ++next;
    } while(curr != circ);

    return next->twin();
  }

  static void approximate_vertex(Execution_context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_at_open_boundary()) {
      return;
    }
    ctx.bounded_approx_pt(vh);

    auto portals_it = ctx->feature_portals.find(vh);
    if(portals_it == ctx->feature_portals.end()) {
      return;
    }
    const auto& portals = portals_it->second;
    if(portals.empty()) {
      return;
    }
    CGAL_assertion_msg(portals.size() == 1, "Vertex should have at most one portal.");
    traverse_portal(ctx, portals.front());
  }

  template <typename DirectionTag>
  static void approximate_halfedge(Execution_context& ctx, const Halfedge_const_handle& he, DirectionTag dir_tag) {
    constexpr bool Is_left_to_right = std::is_same_v<DirectionTag, Left_to_right_tag>;
    using Approx_point_it = std::conditional_t<Is_left_to_right, typename Polyline_geom::const_iterator,
                                               typename Polyline_geom::const_reverse_iterator>;
    using Portals_it = std::conditional_t<Is_left_to_right, typename Portal_vector::const_iterator,
                                          typename Portal_vector::const_reverse_iterator>;

    const Polyline_geom& polyline = he->is_fictitious() ? Polyline_geom{} : ctx.bounded_approx_curve(he);
    auto portal_map_it = ctx->feature_portals.find(he);
    const Portal_vector& portals =
        portal_map_it != ctx->feature_portals.end() ? portal_map_it->second : Portal_vector{};

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

    std::merge(points_iter, points_end, portals_iter, portals_end,
               boost::make_function_output_iterator([&ctx](const Point_or_portal& point_or_portal) {
                 if(auto* portal = std::get_if<Portal>(&point_or_portal)) {
                   traverse_portal(ctx, *portal);
                   return;
                 }
                 *ctx.out_it++ = std::get<Approx_point>(point_or_portal);
               }),
               [](const Portal& portal, const Approx_point& pt) {
                 return Is_left_to_right ? portal.second->point().x() < FT(pt.x())
                                         : portal.second->point().x() > FT(pt.x());
               });
  }

  static void approximate_inner_ccb(Execution_context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_isolated()) {
      approximate_vertex(ctx, vh);
      return;
    }
    approximate_ccb<Inner_ccb_tag>(ctx, get_inner_ccb_first_halfedge(vh));
  }

  static void traverse_portal(Execution_context& ctx, const Portal& portal) { // We have a portal.
    const auto& [source, dest] = portal;
    Approx_point dest_pt = ctx->approx_pt(dest->point());
    Approx_point source_pt = source.has_value() ? source.value() : Approx_point(dest_pt.x(), ctx->ymax());

    *ctx.out_it++ = source_pt;
    approximate_inner_ccb(ctx, dest);
    // We come back here after inner ccb is approximated.
    *ctx.out_it++ = source_pt;

    if(!source.has_value()) {
      ctx.passed_fictitious_edge = true;
    }
  }

  template <typename CcbTag, bool Unbounded = false>
  static void approximate_ccb(Execution_context& ctx, Ccb_halfedge_const_circulator start_circ) {
    constexpr bool Is_outer_ccb = std::is_same_v<CcbTag, Outer_ccb_tag>;
    static_assert(Is_outer_ccb || !Unbounded, "Inner CCBs are impossible to be unbounded.");

    // For unbound ccb, we start on a fictitious edge
    if constexpr(Unbounded) {
      while(!start_circ->is_fictitious()) {
        ++start_circ;
      }
    }

    auto circ = start_circ;
    do {
      if(circ->is_fictitious()) {
        ctx.passed_fictitious_edge = true;
      }

      circ->direction() == ARR_LEFT_TO_RIGHT ? approximate_halfedge(ctx, circ, Left_to_right_tag{})
                                             : approximate_halfedge(ctx, circ, Right_to_left_tag{});
      approximate_vertex(ctx, circ->target());
    } while(++circ != start_circ);

    if constexpr(Unbounded) {
      if(ctx.first_pt.has_value()) {
        ctx.patch(ctx.last_pt.value(), ctx.first_pt.value(), ctx.base_out_it);
      } else {
        ctx.patch(ctx.base_out_it);
      }
    }
  }

public:
  Arr_bounded_approximate_face_2(const Bounded_render_context& ctx,
                                 const Bounded_approximate_point_2& bounded_approx_pt,
                                 const Bounded_approximate_curve_2& bounded_approx_curve)
      : m_ctx(ctx)
      , m_patch(ctx.bbox())
      , m_bounded_approx_pt(bounded_approx_pt)
      , m_bounded_approx_curve(bounded_approx_curve) {}

  const Triangulated_face& operator()(const Face_const_handle& fh) const {
    auto [triangulated_face, inserted] = m_ctx.cache.try_emplace(fh);
    if(!inserted) {
      return triangulated_face;
    }

    if(m_ctx.is_cancelled()) {
      return triangulated_face;
    }

    CGAL_precondition_msg(!fh->is_fictitious(), "Cannot approximate a fictitious face.");

    if(!fh->has_outer_ccb()) {
      // The face is the unbounded face of bounded arrangements, we skip approximating any non degenerate features.
      for(auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
        auto circ = *inner_ccb;
        do {
          if(circ->face() != circ->twin()->face()) {
            // Found non degenerate edge, skip.
            continue;
          }
          m_bounded_approx_curve(circ);
        } while(++circ != *inner_ccb);
      }
      for(auto vh = fh->isolated_vertices_begin(); vh != fh->isolated_vertices_end(); ++vh) {
        m_bounded_approx_pt(vh);
      }

      return triangulated_face;
    }

    auto triangulator = Triangulator(m_ctx);
    auto simplifier = Simplifier(triangulator.insert_iterator(), m_ctx.bbox());
    auto ctx =
        Execution_context(m_ctx, m_patch, simplifier.insert_iterator(), m_bounded_approx_pt, m_bounded_approx_curve);

    if(fh->is_unbounded()) {
      approximate_ccb<Outer_ccb_tag, true>(ctx, fh->outer_ccb());
    } else {
      approximate_ccb<Outer_ccb_tag, false>(ctx, fh->outer_ccb());
    }

    simplifier.dump();

    return triangulated_face = std::move(triangulator);
  }

private:
  const Bounded_render_context& m_ctx;
  const Bounded_approximate_point_2& m_bounded_approx_pt;
  const Bounded_approximate_curve_2& m_bounded_approx_curve;
  const Patch m_patch;
};

} // namespace draw_aos
} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H