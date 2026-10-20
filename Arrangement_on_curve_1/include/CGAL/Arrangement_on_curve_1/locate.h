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
    if (arr.empty()) return Result{std::in_place_index<1>, topo.left_boundary_edge()};

    std::size_t idx = arr.binary_search_vertex(p, cmp);

    // Check if p coincides with the vertex at idx
    if (idx < topo.number_of_vertices() && cmp(p, topo.vertex_point(idx)) == CGAL::EQUAL)
      return Result{std::in_place_index<0>, idx};

    // p lies in an edge: either to the left of vertex idx (idx < n),
    // or to the right of the last vertex (idx == n).
    // In both cases the containing edge is the right edge of vertex idx-1,
    // which is the same as the left edge of vertex idx for idx < n.
    // When idx == 0 the arrangement is non-empty but p < all vertices,
    // so p is in left_boundary_edge() — topo.left_edge(0) gives that.
    if (idx == topo.number_of_vertices()) return Result{std::in_place_index<1>, arr.right_edge(idx - 1)};
    return Result{std::in_place_index<1>, topo.left_edge(idx)};
  }
  // --------------------------------------------------------------------------
  // 2. LINEAR PATH: Sequential / Cyclic Topological Walk (O(n))
  //
  // left_boundary_edge() is the uniform entry point for all topology types:
  //   - open/bounded: the leftmost edge in the linear arrangement.
  //   - identified:   the edge whose left end is at the identification point.
  // has_right_vertex() returns false on the spanning edge (empty arrangement).
  // For identified topologies the walk also stops when the next right vertex
  // is the identification point, because Compare_x_1 must not be called on it.
  // --------------------------------------------------------------------------
  else {
    auto v_pnt_map = arr.vertex_point_map();
    auto curr_e = arr.left_boundary_edge();

    // Predicate: is a point the identification point?
    // Compiled away entirely for non-identified arrangements.
    auto is_on_id = [&](const auto& pt) -> bool {
      if constexpr (Arr::is_identified) return arr.geometry_traits_1().is_on_x_identification_1_object()(pt);
      else return false;
    };

    while (arr.has_right_vertex(curr_e)) {
      auto curr_v = arr.right_vertex(curr_e);
      auto curr_pnt = get(v_pnt_map, curr_v);

      // For identified topologies: the identification vertex terminates the
      // cyclic walk — p lies in the edge that wraps back to the seam.
      if (is_on_id(curr_pnt)) break;

      auto res = cmp(p, curr_pnt);

      if (res == CGAL::EQUAL) return Result{std::in_place_index<0>, curr_v};
      if (res == CGAL::SMALLER) return Result{std::in_place_index<1>, curr_e};

      curr_e = arr.right_edge(curr_v);
    }

    // p lies to the right of all existing interior vertices (linear), or
    // between the last interior vertex and the identification point (cyclic).
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
