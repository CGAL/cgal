// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_LOCATE_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_LOCATE_H

#include <variant>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
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
  using Result = typename Arr::Location_result;

  auto cmp = arr.geometry_traits_1().compare_x_1_object();
  auto v_pnt_map = arr.vertex_point_map();

  // Start with the leftmost unbounded edge (-inf, ...)
  auto curr_e = arr.unbounded_edge();

  // Proceed to the right by traversing the topological graph
  while (arr.has_right_vertex(curr_e)) {
    auto curr_v = arr.right_vertex(curr_e);
    auto res = cmp(p, get(v_pnt_map, curr_v));

    if (res == EQUAL) return Result{curr_v};   // the point coincides with a vertex
    if (res == SMALLER) return Result{curr_e}; // the point lies strictly to the left of the current vertex

    // LARGER: the point lies strictly to the right of the current vertex; Move to the right of this vertex
    curr_e = arr.right_edge(curr_v);
  }

  // The point is to the right of all existing vertices (or the arrangement is empty)
  return Result{curr_e};
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
  using Edge_const_descriptor = typename Arr::Edge_const_descriptor;
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

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
