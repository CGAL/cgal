// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_INSERT_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_INSERT_H

#include <iterator>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/locate.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ==========================================
// INSERT
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Vertex_descriptor
insert(Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr, const typename GeometryTraits::Point_1& p) {
  using Arr = Arrangement_on_curve_1<GeometryTraits, TopologyTraits>;
  using Vertex_descriptor = typename Arr::Vertex_descriptor;
  using Edge_descriptor = typename Arr::Edge_descriptor;

  std::cout << "inserting " << p << std::endl;
  auto loc = locate(arr, p);

  // Point already exists → return existing vertex.
  if (std::holds_alternative<Vertex_descriptor>(loc)) return std::get<Vertex_descriptor>(loc);

  auto e = std::get<Edge_descriptor>(loc);

  bool has_left = arr.has_left_vertex(e);
  bool has_right = arr.has_right_vertex(e);

  // Empty arrangement: single unbounded edge.
  if (! has_left && ! has_right) return arr.insert_empty(p);

  // Leftmost unbounded edge (-inf, v_right): p is to the left of everything.
  if (! has_left) return arr.insert_before(arr.right_vertex(e), p);

  // Rightmost unbounded edge (v_left, +inf): p is to the right of everything.
  if (! has_right) return arr.insert_after(arr.left_vertex(e), p);

  // Bounded interior edge (v_left, v_right): split it.
  return arr.split_edge(e, p);
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
