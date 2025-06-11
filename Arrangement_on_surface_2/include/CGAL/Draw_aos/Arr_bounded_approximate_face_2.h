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
    }
    return;
  }

private:
  const Bbox_2 m_bbox;
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
  using Patch_boundary = internal::Patch_boundary;
  using Triangulated_face = Approx_geom_traits::Triangulated_face;

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
                      const Arr_bounded_approximate_point_2& bounded_approx_pt,
                      const Arr_bounded_approximate_curve_2& bounded_approx_curve,
                      const Patch_boundary& patch_boundary)
        : Arr_context_delegator(ctx)
        , triangulator(triangulator)
        , patch_boundary(patch_boundary)
        , bounded_approx_pt(bounded_approx_pt)
        , bounded_approx_curve(bounded_approx_curve) {}

  public:
    const Arr_bounded_approximate_point_2& bounded_approx_pt;
    const Arr_bounded_approximate_curve_2& bounded_approx_curve;
    const Patch_boundary& patch_boundary;
    Arr_bounded_face_triangulator& triangulator;
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

  template <typename Ccb_tag>
  static void approximate_ccb(Execution_context& ctx, const Ccb_halfedge_const_circulator& start_circ) {
    constexpr bool Is_outer_ccb = std::is_same_v<Ccb_tag, Outer_ccb_tag>;
    auto ccb_constraint = ctx.triangulator.make_ccb_constraint<Ccb_tag>();
    auto out_it = ccb_constraint.insert_iterator();

    std::optional<Approx_point> ccb_last_pt, ccb_first_pt;
    auto counter_clockwise_start_circ = Is_outer_ccb ? start_circ : Ccb_halfedge_const_circulator(start_circ->twin());
    auto circ = counter_clockwise_start_circ;
    do {
      bool is_he_first_pt = true;
      auto patch_out_it = boost::make_function_output_iterator([&](const Approx_point& pt) {
        if(ccb_last_pt == pt || !ctx->contains(pt)) {
          return;
        }
        if(is_he_first_pt && ccb_last_pt.has_value()) {
          ctx.patch_boundary(ccb_last_pt.value(), pt, out_it);
        }

        *out_it++ = pt;
        ccb_last_pt = pt;
        if(!ccb_first_pt.has_value()) {
          ccb_first_pt = pt;
        }
        is_he_first_pt = false;
      });

      approximate_halfedge_of_ccb(ctx, circ, patch_out_it);
      approximate_vertex(ctx, circ->target());
    } while(++circ != counter_clockwise_start_circ);

    if(Is_outer_ccb && !ccb_first_pt.has_value()) {
      *out_it++ = Approx_point(ctx->xmin(), ctx->ymin());
      *out_it++ = Approx_point(ctx->xmax(), ctx->ymin());
      *out_it++ = Approx_point(ctx->xmax(), ctx->ymax());
      *out_it++ = Approx_point(ctx->xmin(), ctx->ymax());
    }
    if(ccb_first_pt.has_value() && ccb_first_pt != ccb_last_pt) {
      // Close the ccb
      ctx.patch_boundary(ccb_last_pt.value(), ccb_first_pt.value(), out_it);
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
      // The face is the unbounded face of bounded arrangements
      return triangulated_face;
    }

    Arr_bounded_face_triangulator triangulator(m_ctx);
    Execution_context ctx(m_ctx, triangulator, m_point_approx, m_curve_approx, m_patch_boundary);

    approximate_ccb<Outer_ccb_tag>(ctx, fh->outer_ccb());
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