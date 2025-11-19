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

#ifndef CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H
#define CGAL_DRAW_AOS_ARR_APPROXIMATION_CACHE_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <boost/range/iterator_range.hpp>

#include <CGAL/Arr_enums.h>
#include <CGAL/unordered_flat_map.h>
#include <CGAL/Draw_aos/type_utils.h>

namespace CGAL {
namespace draw_aos {

/** @brief Cache class for approximating arrangement on surface.
 *
 * When iterating over the arrangement dcel, a feature(vertex, halfedge, face) might be visited multiple times.
 * This cache stores the approximated geometry for each feature to avoid redundant calculations.
 * @tparam Arrangement
 */
template <typename Arrangement>
class Arr_approximation_cache {
  using Geom_traits = typename Arrangement::Geometry_traits_2;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;

  using Vertex_cache_obj = typename Approx_traits::Point;
  using Halfedge_cache_obj = typename Approx_traits::Polyline;
  using Face_cache_obj = typename Approx_traits::Triangle_soup;

  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_iterator;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_cache = unordered_flat_map<Vertex_const_handle, Vertex_cache_obj>;
  using Halfedge_cache = unordered_flat_map<Halfedge_const_handle, Halfedge_cache_obj>;
  using Face_cache = unordered_flat_map<Face_const_handle, Face_cache_obj>;

public:
  Arr_approximation_cache() = default;

  const Vertex_cache& vertices() const { return m_vertices; }
  const Halfedge_cache& halfedges() const { return m_halfedges; }
  const Face_cache& faces() const { return m_faces; }

  Vertex_cache& vertices() { return m_vertices; }
  Halfedge_cache& halfedges() { return m_halfedges; }
  Face_cache& faces() { return m_faces; }

private:
  Vertex_cache m_vertices;
  Halfedge_cache m_halfedges;
  Face_cache m_faces;
};

} // namespace draw_aos
} // namespace CGAL

#endif
