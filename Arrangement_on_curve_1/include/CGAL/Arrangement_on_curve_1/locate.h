// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
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
// Uses \f$O(log n)\f$O binary search when `BinarySearch` = `true`,
// otherwise falls back to \f$O(n)\f$O topological graph traversal.
// ==========================================
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>::Location_result
locate_impl(Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr,
            const typename GeometryTraits_1::Point_1& p) {
  using Aoc = Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>;
  using Result = typename Aoc::Location_result;

  auto cmp = arr.geometry_traits_1().compare_x_1_object();

  // --------------------------------------------------------------------------
  // 1. FAST PATH: Binary Search (\f$O(log n)\f$O)
  // --------------------------------------------------------------------------
  if constexpr (Aoc::binary_search_enabled) {
    auto& topo = arr.topology_traits();
    if (arr.empty()) return Result{std::in_place_index<1>, topo.left_boundary_edge()};

    std::size_t idx = arr.binary_search_vertex(p, cmp);

    // Check if p coincides with the vertex at idx
    if (idx < topo.number_of_vertices() && cmp(p, topo.vertex_point(idx)) == EQUAL)
      return Result{std::in_place_index<0>, idx};

    // p lies in an edge: either to the left of vertex idx (idx < n),
    // or to the right of the last vertex (idx == n).
    // In both cases the containing edge is the right edge of vertex idx-1,
    // which is the same as the left edge of vertex idx for idx < n.
    // When idx == 0 the arrangement is non-empty but p < all vertices,
    // so p is in left_boundary_edge(); topo.left_edge(0) gives that.
    if (idx == topo.number_of_vertices()) return Result{std::in_place_index<1>, arr.right_edge(idx - 1)};
    return Result{std::in_place_index<1>, topo.left_edge(idx)};
  }
  // --------------------------------------------------------------------------
  // 2. LINEAR PATH: Sequential / Cyclic Topological Walk (\f$O(n)\f$O)
  //
  // left_boundary_edge() is the uniform entry point for all topology types:
  //   - open/closed: the leftmost edge in the linear arrangement.
  //   - identified: the edge whose left end is the identification vertex, or,
  //     if it does not exist, the edge that wraps around the identification
  //     point.
  // has_right_vertex() returns false on the spanning edge (empty arrangement).
  // --------------------------------------------------------------------------
  else {
    auto v_pnt_map = arr.vertex_point_map();
    auto curr_e = arr.left_boundary_edge();

    if constexpr (Aoc::is_identified) {
      // The identification point is handled first, because Compare_x_1 must
      // not be called on it.  If the identification vertex does not exist,
      // the point lies in left_boundary_edge().
      auto& topo = arr.topology_traits();
      const bool has_id_vertex = topo.has_identification_vertex();
      if (arr.geometry_traits_1().is_on_x_identification_1_object()(p)) {
        if (has_id_vertex) return Result{std::in_place_index<0>, topo.identification_vertex()};
        return Result{std::in_place_index<1>, curr_e};
      }

      // Every edge has a right vertex, so the walk is bounded by the number of
      // the remaining vertices, which it visits in increasing order.  It ends
      // at the edge that precedes the identification vertex, or, if it does
      // not exist, back at left_boundary_edge().
      for (auto n = arr.number_of_vertices() - (has_id_vertex ? 1 : 0); n > 0; --n) {
        auto curr_v = arr.right_vertex(curr_e);
        auto res = cmp(p, get(v_pnt_map, curr_v));
        if (res == EQUAL) return Result{std::in_place_index<0>, curr_v};
        if (res == SMALLER) return Result{std::in_place_index<1>, curr_e};
        curr_e = arr.right_edge(curr_v);
      }

      // p lies between the last vertex and the identification point.
      return Result{std::in_place_index<1>, curr_e};
    }
    else {
      // For a closed left side: the left boundary vertex if exists is the left vertex of the left boundary edge.
      if constexpr (Aoc::left_is_closed) {
        if (arr.has_left_vertex(curr_e)) {
          auto v = arr.left_vertex(curr_e);
          if (cmp(p, get(v_pnt_map, v)) == EQUAL) return Result{std::in_place_index<0>, v};
        }
      }

      while (arr.has_right_vertex(curr_e)) {
        auto curr_v = arr.right_vertex(curr_e);
        auto res = cmp(p, get(v_pnt_map, curr_v));
        if (res == EQUAL) return Result{std::in_place_index<0>, curr_v};
        if (res == SMALLER) return Result{std::in_place_index<1>, curr_e};
        curr_e = arr.right_edge(curr_v);
      }

      // p lies to the right of all existing vertices.
      return Result{std::in_place_index<1>, curr_e};
    }
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

  using Aoc = Arrangement_on_curve_1<Geometry_traits_1, Topology_traits, BinarySearch>;
  using Vertex_const_descriptor = typename Aoc::Vertex_const_descriptor;
  using Edge_const_descriptor = typename Aoc::Edge_const_descriptor;
  using Result = typename Aoc::Const_location_result;

  auto mutable_result = locate_impl(const_cast<Aoc&>(arr), p);

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
