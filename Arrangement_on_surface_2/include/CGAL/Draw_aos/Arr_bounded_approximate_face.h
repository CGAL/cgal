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

#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_FACE_H
#include <iterator>
#include <optional>
#include <utility>
#include <algorithm>

#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_halfedge.h>
#include <CGAL/Draw_aos/Arr_bounded_approximate_vertex.h>
#include <CGAL/Draw_aos/Arr_bounded_face_triangulator.h>
#include <CGAL/Draw_aos/Arr_render_context.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {

namespace draw_aos {

/*! \brief Functor to approximate arrangement face with triangles within a bounding box.
 *
 * \tparam Arrangement
 */
template <typename Arrangement>
class Arr_bounded_approximate_face {
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Ccb_halfedge_const_circulator = typename Arrangement::Ccb_halfedge_const_circulator;

  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Point = typename Approx_traits::Point;
  using Polyline = typename Approx_traits::Polyline;
  using Triangle_soup = typename Approx_traits::Triangle_soup;

  using Bounded_approximate_vertex = Arr_bounded_approximate_vertex<Arrangement>;
  using Bounded_approximate_halfedge = Arr_bounded_approximate_halfedge<Arrangement>;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;
  using Triangulator = Arr_bounded_face_triangulator<Arrangement>;

  static constexpr bool Is_on_curved_surface = is_or_derived_from_curved_surf_traits_v<Geom_traits>;

  struct Left_to_right_tag {};
  struct Right_to_left_tag {};

private:
  /*! \brief A stateful geometry simplifier that simplifies horizontal and vertical segments
   *
   * \tparam OutputIterator
   */
  template <typename OutputIterator>
  class Colinear_simplifier {
  public:
    Colinear_simplifier(OutputIterator out_it) : m_out_it(out_it) {}

    void dump() {
      if (m_start.has_value()) {
        *m_out_it++ = m_start.value();
        m_start.reset();
      }
      if (m_mid.has_value()) {
        *m_out_it++ = m_mid.value();
        m_mid.reset();
      }
    }

    void push_back(Point p) {
      if (m_mid.has_value()) {
        if (((p.y() == m_mid->y()) && (p.y() == m_start->y())) || ((p.x() == m_mid->x()) && (p.x() == m_start->x())))
          // Three points are collinear horizontally or vertically.
          m_mid = p;
        else {
          *m_out_it++ = m_start.value();
          m_start = m_mid;
          m_mid = p;
        }
        return;
      }

      if (m_start.has_value())
        m_mid = p;
      else
        m_start = p;
    }

    ~Colinear_simplifier() { dump(); }

  private:
    OutputIterator m_out_it;
    std::optional<Point> m_start, m_mid;
  };

  class Context : public Bounded_render_context {
    using Simplifier = Colinear_simplifier<std::back_insert_iterator<Triangulator>>;

  public:
    Context(const Bounded_render_context& ctx, Triangulator& triangulator) :
      Bounded_render_context(ctx),
      m_triangulator(triangulator) {
      if constexpr(!Is_on_curved_surface) m_simplifier.emplace(std::back_inserter(m_triangulator));
    }

    // Let's not accidentally copy this object.
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

    void insert(Point pt) {
      if (Approx_traits::is_null(pt) || pt == m_last_pt) return;
      pt = Point(pt.x(), std::clamp(pt.y(), this->ymin(), this->ymax()));
      if constexpr(!Is_on_curved_surface) {
        m_simplifier->push_back(pt);
        return;
      }
      m_triangulator.push_back(pt);
      m_last_pt = pt;
    }

    void start_ccb() { m_triangulator.start_constraint(); }

    void end_ccb() {
      if constexpr(!Is_on_curved_surface) m_simplifier->dump();
      m_triangulator.end_constraint();
    }

    const std::optional<Point>& last_pt() const { return m_last_pt; }

  private:
    Triangulator& m_triangulator;
    // Colinear simplifier is only used for optimizing planar arrangements.
    std::optional<Simplifier> m_simplifier;
    std::optional<Point> m_last_pt;
  };

private:
  static Arr_parameter_space side_of_fict_edge(const Halfedge_const_handle& he) {
    const auto& source = he->source();
    const auto& target = he->target();
    auto sx = source->parameter_space_in_x();
    auto sy = source->parameter_space_in_y();
    auto tx = target->parameter_space_in_x();
    auto ty = target->parameter_space_in_y();
    if (sx == tx && sx != ARR_INTERIOR) return sx;
    if (sy == ty && sy != ARR_INTERIOR) return sy;
    CGAL_assertion(false && "Unexpected parameter space for fictitious edge ends.");
    return ARR_INTERIOR;
  }

  // Generate dummy segment for fictitious edge he at its corresponding boundary.
  static Polyline approximate_fict_edge(const Context& ctx, const Halfedge_const_handle& he) {
    auto side = side_of_fict_edge(he);
    // There's no need to handle fictitious edges on left or right boundaries.
    if (side == ARR_LEFT_BOUNDARY || side == ARR_RIGHT_BOUNDARY) return Polyline{};
    if (side == ARR_BOTTOM_BOUNDARY) return Polyline{ctx.bottom_left(), ctx.bottom_right()};
    if (side == ARR_TOP_BOUNDARY) return Polyline{ctx.top_right(), ctx.top_left()};
    CGAL_assertion(false && "Unexpected side for a fictitious edge.");
    return Polyline{};
  }

  void approximate_vertex(Context& ctx, const Vertex_const_handle& vh) const {
    if (vh->is_at_open_boundary()) return;
    m_bounded_approx_vertex(vh);
  }

  void approximate_halfedge(Context& ctx, const Halfedge_const_handle& he) const {
    const Polyline& polyline = he->is_fictitious() ? approximate_fict_edge(ctx, he) : m_bounded_approx_halfedge(he);
    for (const auto& curr_pt : polyline) ctx.insert(curr_pt);
  }

  void approximate_ccb(Context& ctx, Ccb_halfedge_const_circulator start) const {
    // Try to start on a concrete halfedge.
    // For any unbounded face, there can't be more than 4 adjacent fictitious edges.
    for (int i = 0; i < 4 && start->is_fictitious(); ++i) ++start;

    ctx.start_ccb();
    auto circ = start;
    do {
      approximate_halfedge(ctx, circ);
      approximate_vertex(ctx, circ->target());
    } while(++circ != start);
    ctx.end_ccb();
  }

public:
  Arr_bounded_approximate_face(const Bounded_render_context& ctx) :
    m_ctx(ctx),
    m_bounded_approx_halfedge(ctx),
    m_bounded_approx_vertex(ctx)
  {}

  /*! \brief Approximate an arrangement face with a bunch of triangles.
   *
   * \param fh
   * \return const Triangulated_face&
   */
  const Triangle_soup& operator()(const Face_const_handle& fh) const {
    CGAL_precondition_msg(!fh->is_fictitious(), "Cannot approximate a fictitious face.");

    auto [iter, inserted] = m_ctx.m_cache.faces().try_emplace(fh);
    Triangle_soup& ts = iter->second;
    if (! inserted || m_ctx.is_cancelled()) return ts;
    auto triangulator = Triangulator(m_ctx, fh);
    auto ctx = Context(m_ctx, triangulator);

    if (! Is_on_curved_surface && !fh->has_outer_ccb()) {
      // Skip approximation of the unbounded face in planar arrangements.
      // However, degenerate holes still need to be approximated.
      for (auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb) {
        auto circ = *inner_ccb;
        do {
          if (circ->face() != circ->twin()->face()) continue;
          m_bounded_approx_halfedge(circ);
        } while(++circ != *inner_ccb);
      }
      for (auto vh = fh->isolated_vertices_begin(); vh != fh->isolated_vertices_end(); ++vh) m_bounded_approx_vertex(vh);
      return ts;
    }

    for (auto outer_ccb = fh->outer_ccbs_begin(); outer_ccb != fh->outer_ccbs_end(); ++outer_ccb)
      approximate_ccb(ctx, *outer_ccb);
    for (auto inner_ccb = fh->inner_ccbs_begin(); inner_ccb != fh->inner_ccbs_end(); ++inner_ccb)
      approximate_ccb(ctx, *inner_ccb);
    for (auto iso_vertex = fh->isolated_vertices_begin(); iso_vertex != fh->isolated_vertices_end(); ++iso_vertex)
      approximate_vertex(ctx, iso_vertex);

    return ts = std::move(triangulator);
  }

private:
  const Bounded_render_context& m_ctx;
  const Bounded_approximate_halfedge m_bounded_approx_halfedge;
  const Bounded_approximate_vertex m_bounded_approx_vertex;
};

} // namespace draw_aos
} // namespace CGAL

#endif
