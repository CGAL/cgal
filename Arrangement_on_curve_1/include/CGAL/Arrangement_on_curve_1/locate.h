// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
// Uses O(log n) binary search when BinarySearch = true,
// otherwise falls back to O(n) topological graph traversal.
// ==========================================
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>::Location_result
locate_impl(Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr,
            const typename GeometryTraits_1::Point_1& p) {
  using Arr = Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>;
  using Result = typename Arr::Location_result;

  auto cmp = arr.geometry_traits_1().compare_x_1_object();

  // --------------------------------------------------------------------------
  // 1. FAST PATH: Binary Search (O(log n))
  // --------------------------------------------------------------------------
  if constexpr (Arr::binary_search_enabled) {
    auto& topo = arr.topology_traits();
    if (arr.empty()) return Result{std::in_place_index<1>, topo.unbounded_left_edge()};

    std::size_t idx = arr.binary_search_vertex(p, cmp);

    // p is to the right of all existing vertices
    if (idx == topo.number_of_vertices()) return Result{std::in_place_index<1>, topo.unbounded_right_edge()};

    // Check if p coincides with the vertex at idx
    if (cmp(p, topo.vertex_point(idx)) == CGAL::EQUAL) return Result{std::in_place_index<0>, idx};

    // p is strictly to the left of vertex at idx;
    // thus, p lies in the edge immediately to the left of vertex idx.
    return Result{std::in_place_index<1>, topo.left_edge(idx)};
  }
  // --------------------------------------------------------------------------
  // 2. LINEAR PATH: Sequential Topological Walk (O(n))
  // --------------------------------------------------------------------------
  else {
    auto v_pnt_map = arr.vertex_point_map();

    // Start with the leftmost unbounded edge (-inf, ...)
    auto curr_e = arr.unbounded_left_edge();

    // Proceed to the right by traversing the topological graph
    while (arr.has_right_vertex(curr_e)) {
      auto curr_v = arr.right_vertex(curr_e);
      auto res = cmp(p, get(v_pnt_map, curr_v));

      if (res == CGAL::EQUAL) return Result{std::in_place_index<0>, curr_v}; // point coincides with vertex
      if (res == CGAL::SMALLER) return Result{std::in_place_index<1>, curr_e}; // point lies strictly to left of current vertex

      // LARGER: move to right edge of this vertex
      curr_e = arr.right_edge(curr_v);
    }

    // The point is to the right of all existing vertices (or arrangement is empty)
    return Result{std::in_place_index<1>, curr_e};
  }
}

// ==========================================
// LOCATE (Mutable)
// ==========================================
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>::Location_result
locate(Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits_1::Point_1& p)
{ return locate_impl(arr, p); }

// ==========================================
// LOCATE (Const)
// Delegates to mutable impl via a const_cast (safe: does not mutate).
// ==========================================
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>::Const_location_result
locate(const Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits_1::Point_1& p) {
  using Geometry_traits_1 = GeometryTraits_1;
  using Topology_traits = TopologyTraits;

  using Arr = Arrangement_on_curve_1<Geometry_traits_1, Topology_traits, BinarySearch>;
  using Vertex_const_descriptor = typename Arr::Vertex_const_descriptor;
  using Edge_const_descriptor = typename Arr::Edge_const_descriptor;
  using Result = typename Arr::Const_location_result;

  auto mutable_result = locate_impl(const_cast<Arr&>(arr), p);

  if (mutable_result.index() == 0) {
    auto v = std::get<0>(mutable_result);
    return Result{std::in_place_index<0>, Vertex_const_descriptor{v}};
  }

  auto e = std::get<1>(mutable_result);
  return Result{std::in_place_index<1>, Edge_const_descriptor{e}};
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
