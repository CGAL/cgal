// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_FUNCTIONS_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_FUNCTIONS_H

#include <variant>
#include <iterator>

#include <CGAL/Arrangement_on_curve_1.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ==========================================
// Internal helper: locate, returning mutable descriptors.
// Walks the sorted vertex sequence and finds where p falls.
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Location_result
locate_impl(Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr, const typename GeometryTraits::Point_1& p) {
  using Arr = Arrangement_on_curve_1<GeometryTraits, TopologyTraits>;
  using Vertex_descriptor = typename Arr::Vertex_descriptor;
  using Edge_descriptor = typename Arr::Edge_descriptor;
  using Result = typename Arr::Location_result;

  auto cmp = arr.geometry_traits_1().compare_x_1_object();
  auto v_pnt_map = arr.vertex_point_map();
  auto& topo = arr.topology_traits();

  // Walk through the mutable vertex list.
  for (auto vit = topo.vertices_begin(); vit != topo.vertices_end(); ++vit) {
    auto res = cmp(p, get(v_pnt_map, vit));

    if (res == EQUAL) return Result{Vertex_descriptor{vit}};
    if (res == SMALLER) return Result{Edge_descriptor{arr.left_edge(vit)}};
    // LARGER: continue walking right
  }

  // p is beyond all vertices → rightmost edge (or the sole unbounded edge).
  if (topo.is_empty()) return Result{Edge_descriptor{arr.unbounded_edge()}};

  return Result{Edge_descriptor{topo.right_edge(std::prev(topo.vertices_end()))}};
}

// ==========================================
// LOCATE (Mutable)
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Location_result
locate(Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr, const typename GeometryTraits::Point_1& p)
{ return locate_impl(arr, p); }

// ==========================================
// LOCATE (Const)
// Delegates to the mutable impl via a const_cast (safe: we do not mutate).
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Const_location_result
locate(const Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr, const typename GeometryTraits::Point_1& p) {
  using Arr = Arrangement_on_curve_1<GeometryTraits, TopologyTraits>;
  using Vertex_const_descriptor = typename Arr::Vertex_const_descriptor;
  using Edge_const_descriptor   = typename Arr::Edge_const_descriptor;
  using Result = typename Arr::Const_location_result;

  auto mutable_result = locate_impl(const_cast<Arrangement_on_curve_1<GeometryTraits, TopologyTraits>&>(arr), p);

  // Convert mutable descriptors → const descriptors.
  using Vertex_descriptor = typename Arr::Vertex_descriptor;
  using Edge_descriptor = typename Arr::Edge_descriptor;

  if (std::holds_alternative<Vertex_descriptor>(mutable_result))
    return Result{Vertex_const_descriptor{std::get<Vertex_descriptor>(mutable_result)}};
  if (std::holds_alternative<Edge_descriptor>(mutable_result))
    return Result{Edge_const_descriptor{std::get<Edge_descriptor>(mutable_result)}};
  return Result{static_cast<void*>(nullptr)};
}

// ==========================================
// INSERT
// ==========================================
template <typename GeometryTraits, typename TopologyTraits>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits>::Vertex_descriptor
insert(Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr, const typename GeometryTraits::Point_1& p) {
  using Arr = Arrangement_on_curve_1<GeometryTraits, TopologyTraits>;
  using Vertex_descriptor = typename Arr::Vertex_descriptor;
  using Edge_descriptor = typename Arr::Edge_descriptor;

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
