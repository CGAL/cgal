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

#ifndef CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_VERTEX_H
#define CGAL_DRAW_AOS_ARR_BOUNDED_APPROXIMATE_VERTEX_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/Draw_aos/type_utils.h>
#include <CGAL/Draw_aos/Arr_render_context.h>

namespace CGAL {
namespace draw_aos {

template <typename Arrangement>
class Arr_bounded_approximate_vertex {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Point_2 = typename Geom_traits::Point_2;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Point_geom = typename Arr_approximate_traits<Geom_traits>::Point;
  using Bounded_render_context = Arr_bounded_render_context<Arrangement>;

public:
  Arr_bounded_approximate_vertex(const Bounded_render_context& ctx) : m_ctx(ctx) {}

  /** @brief Approximate a vertex within the x-bounded range.
   *
   * The function uses cached values if available.
   * @precondition: The vertex must have an associated point.
   *
   * @param vh the vertex handle
   * @return const Point_geom&
   */
  const Point_geom& operator()(const Vertex_const_handle& vh) const {
    auto [iter, inserted] = m_ctx.m_cache.vertices().try_emplace(vh);
    Point_geom& point = iter->second;
    if (! inserted) return point;
    return point = m_ctx.to_uv(m_ctx.m_traits.approximate_2_object()(vh->point()));
  }

private:
  const Bounded_render_context& m_ctx;
};

} // namespace draw_aos
} // namespace CGAL

#endif
