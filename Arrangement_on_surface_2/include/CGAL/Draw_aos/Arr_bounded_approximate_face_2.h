#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H

#include "CGAL/Arr_enums.h"
#include "CGAL/Bbox_2.h"
#include "CGAL/Draw_aos/Arr_bounded_approximate_curve_2.h"
#include "CGAL/Draw_aos/Arr_bounded_approximate_point_2.h"
#include "CGAL/Draw_aos/Arr_bounded_face_triangulator.h"
#include "CGAL/Draw_aos/Arr_render_context.h"
#include "CGAL/basic.h"
#include <CGAL/Draw_aos/Arr_approximation_geometry_traits.h>
#include <CGAL/Draw_aos/helpers.h>
#include <algorithm>
#include <boost/iterator/function_output_iterator.hpp>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <optional>
#include <type_traits>

namespace CGAL {

namespace internal {
/**
 * @brief Patches corners between two boundary points of the bbox
 * counter-clockwisely.
 */
class Patch_boundary
{
  using Approx_point = Arr_approximation_geometry_traits::Approx_point;

  enum class Side_of_boundary {
    Top = 0,
    Left = 1,
    Bottom = 2,
    Right = 3,
    None = -1,
  };

private:
  Side_of_boundary side_of_boundary(const Approx_point& pt) const {
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

  bool is_on_boundary(const Approx_point& pt) const { return side_of_boundary(pt) != Side_of_boundary::None; }

  double distance_on_same_side(const Approx_point& pt1, const Approx_point& pt2) const {
    return std::abs(pt1.x() - pt2.x()) + std::abs(pt1.y() - pt2.y());
  }

public:
  Patch_boundary(Bbox_2 bbox)
      : m_bbox(bbox) {}

  template <typename OutputIterator>
  void operator()(const Approx_point& from, const Approx_point& to, OutputIterator out_it) const {
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
      std::cout << "Patching corner: " << corner << std::endl;
    }
    return;
  }

private:
  const Bbox_2 m_bbox;
};

template <typename OutputIterator>
class Geom_simplifier
{
  using Approx_point = Arr_approximation_geometry_traits::Approx_point;

private:
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

public:
  Geom_simplifier(OutputIterator& out_it, const Bbox_2& bbox)
      : m_out_it(out_it)
      , m_bbox(bbox) {}

  decltype(auto) insert_iterator() {
    return boost::make_function_output_iterator([this](const Approx_point& p) {
      if(m_mid.has_value()) {
        if(m_mid.value() == p) {
          return;
        }
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
        if(m_start.value() == p) {
          return;
        }
        m_mid = p;
      } else {
        m_start = p;
      }
    });
  }

  ~Geom_simplifier() { dump(); }

private:
  OutputIterator& m_out_it;
  std::optional<Approx_point> m_start, m_mid;
  Bbox_2 m_bbox;
};

} // namespace internal

/**
 * @brief Bounded face approximation for arrangements.
 * @note Member functions are not thread-safe.
 */
class Arr_bounded_approximate_face_2
{
  using Approx_geom_traits = Arr_approximation_geometry_traits;
  using Face_const_handle = Arrangement::Face_const_handle;
  using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
  using Vertex_const_handle = Arrangement::Vertex_const_handle;
  using Polyline_geom = Approx_geom_traits::Polyline_geom;
  using Ccb_halfedge_const_circulator = Arrangement::Ccb_halfedge_const_circulator;
  using Approx_point = Approx_geom_traits::Approx_point;
  using Triangulated_face = Approx_geom_traits::Triangulated_face;
  using Patch_boundary = internal::Patch_boundary;

  template <typename OutputIterator>
  using Geom_simplifier = internal::Geom_simplifier<OutputIterator>;

  struct Left_to_right_tag
  {};
  struct Right_to_left_tag
  {};

private:
  class Execution_context : public Arr_context_delegator<Arr_bounded_render_context>
  {
  public:
    Execution_context(const Arr_bounded_render_context& ctx,
                      Arr_bounded_face_triangulator& triangulator,
                      const Patch_boundary& patch_boundary,
                      const Arr_bounded_approximate_point_2& bounded_approx_pt,
                      const Arr_bounded_approximate_curve_2& bounded_approx_curve)
        : Arr_context_delegator(ctx)
        , triangulator(triangulator)
        , patch_boundary(patch_boundary)
        , bounded_approx_pt(bounded_approx_pt)
        , bounded_approx_curve(bounded_approx_curve) {}

  public:
    const Arr_bounded_approximate_point_2& bounded_approx_pt;
    const Arr_bounded_approximate_curve_2& bounded_approx_curve;
    Arr_bounded_face_triangulator& triangulator;
    const Patch_boundary& patch_boundary;
  };

private:
  static void approximate_vertex(Execution_context& ctx, const Vertex_const_handle& vh) {
    if(vh->is_at_open_boundary()) {
      return;
    }
    ctx.bounded_approx_pt(vh);
  }

  template <typename OutputIterator>
  static void
  approximate_halfedge_of_ccb(Execution_context& ctx, const Halfedge_const_handle& he, OutputIterator& out_it) {
    if(he->is_fictitious()) {
      return;
    }

    const Polyline_geom& polyline = ctx.bounded_approx_curve(he);
    if(he->direction() == ARR_LEFT_TO_RIGHT) {
      std::copy(polyline.begin(), polyline.end(), out_it);
    } else {
      std::copy(polyline.rbegin(), polyline.rend(), out_it);
    }
  }

  template <typename CcbTag, bool Bounded = true>
  static void approximate_ccb(Execution_context& ctx, Ccb_halfedge_const_circulator start_circ) {
    constexpr bool Is_outer_ccb = std::is_same_v<CcbTag, Outer_ccb_tag>;
    static_assert(Is_outer_ccb || Bounded, "Inner CCBs are impossible to be unbounded.");

    // For unbound ccb, we start on a fictitious edge
    if constexpr(!Bounded) {
      while(!start_circ->is_fictitious()) {
        ++start_circ;
      }
    }

    auto ccb_constraint = ctx.triangulator.make_ccb_constraint<CcbTag>();
    auto constraint_out_it = ccb_constraint.insert_iterator();

    auto simplifier = Geom_simplifier(constraint_out_it, ctx->bbox());
    auto simplifier_out_it = simplifier.insert_iterator();

    auto circ = start_circ;
    std::optional<Approx_point> last_pt;

    // These vars are used only in unbounded ccb.
    std::optional<Approx_point> first_pt;
    bool passed_fictitious_edge = false;

    auto he_process_out_it = boost::make_function_output_iterator([&](const Approx_point& pt) {
      Approx_point regulated_pt(pt.x(), std::clamp(pt.y(), ctx->ymin(), ctx->ymax()));
      if(last_pt == regulated_pt) {
        return;
      }

      *simplifier_out_it++ = regulated_pt;

      if constexpr(!Bounded) {
        // TODO: nesting too deep and looks ugly
        if(passed_fictitious_edge) {
          passed_fictitious_edge = false;
          if(last_pt.has_value()) {
            ctx.patch_boundary(last_pt.value(), regulated_pt, simplifier_out_it);
          }
        }
        if(!first_pt.has_value()) {
          first_pt = regulated_pt;
        }
      }

      last_pt = regulated_pt;
    });

    do {
      if constexpr(!Bounded) {
        if(circ->is_fictitious()) {
          passed_fictitious_edge = true;
        }
      }
      approximate_halfedge_of_ccb(ctx, circ, he_process_out_it);
      approximate_vertex(ctx, circ->target());
    } while(++circ != start_circ);

    if constexpr(!Bounded) {
      if(!first_pt.has_value()) {
        *simplifier_out_it++ = Approx_point(ctx->xmin(), ctx->ymin());
        *simplifier_out_it++ = Approx_point(ctx->xmin(), ctx->ymax());
        *simplifier_out_it++ = Approx_point(ctx->xmax(), ctx->ymax());
        *simplifier_out_it++ = Approx_point(ctx->xmax(), ctx->ymin());

      } else {
        ctx.patch_boundary(last_pt.value(), first_pt.value(), simplifier_out_it);
      }
    }
  }

public:
  Arr_bounded_approximate_face_2(const Arr_bounded_render_context& ctx,
                                 const Arr_bounded_approximate_point_2& point_approx,
                                 const Arr_bounded_approximate_curve_2& curve_approx)
      : m_ctx(ctx)
      , m_patch_boundary(ctx.bbox())
      , m_point_approx(point_approx)
      , m_curve_approx(curve_approx) {}

  const Triangulated_face& operator()(const Face_const_handle& fh) const {
    auto [triangulated_face, inserted] = m_ctx.cache.try_emplace(fh);
    if(!inserted) {
      return triangulated_face;
    }

    if(m_ctx.is_cancelled()) {
      return triangulated_face;
    }

    CGAL_assertion_msg(!fh->is_fictitious(), "Cannot approximate a fictitious face.");

    if(!fh->has_outer_ccb()) {
      // The face is the unbounded face of bounded arrangements, we skip approximating any non degenerate features.
      for(auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
        if((*inner_ccb)->twin()->face() != (*inner_ccb)->face()) {
          continue;
        }
        m_curve_approx(*inner_ccb);
      }
      for(auto isolated_vh = fh->isolated_vertices_begin(); isolated_vh != fh->isolated_vertices_end(); ++isolated_vh) {
        m_point_approx(isolated_vh);
      }
      return triangulated_face;
    }

    Arr_bounded_face_triangulator triangulator(m_ctx);
    Execution_context ctx(m_ctx, triangulator, m_patch_boundary, m_point_approx, m_curve_approx);

    if(fh->is_unbounded()) {
      approximate_ccb<Outer_ccb_tag, false>(ctx, fh->outer_ccb());
    } else {
      approximate_ccb<Outer_ccb_tag, true>(ctx, fh->outer_ccb());
    }

    for(auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
      approximate_ccb<Inner_ccb_tag>(ctx, *inner_ccb);
    }

    for(auto isolated_vh = fh->isolated_vertices_begin(); isolated_vh != fh->isolated_vertices_begin(); ++isolated_vh) {
      approximate_vertex(ctx, isolated_vh);
    }

    return triangulated_face = (std::move(triangulator));
  }

private:
  const Arr_bounded_render_context& m_ctx;
  const Arr_bounded_approximate_point_2& m_point_approx;
  const Arr_bounded_approximate_curve_2& m_curve_approx;
  const Patch_boundary m_patch_boundary;
};
} // namespace CGAL
#endif // CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_2_H