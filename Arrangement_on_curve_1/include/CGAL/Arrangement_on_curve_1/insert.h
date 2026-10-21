// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_INSERT_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_INSERT_H

#include <iterator>
#include <variant>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/locate.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! inserts a point into an arrangement on a curve.
 */
template <typename GeometryTraits, typename TopologyTraits, bool BinarySearch>
typename Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>::Vertex_descriptor
insert(Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>& arr,
       const typename GeometryTraits::Point_1& p) {
  using Aoc = Arrangement_on_curve_1<GeometryTraits, TopologyTraits, BinarySearch>;
  using Vertex_descriptor = typename Aoc::Vertex_descriptor;
  using Edge_descriptor = typename Aoc::Edge_descriptor;

  auto loc = locate(arr, p);

  // 1. In Vector Mode (`UseVector` = true): Location_result is `std::variant<std::size_t, std::size_t>`.
  // Checking `loc.index() == 0` and retrieving `std::get<0>(loc)` or `std::get<1>(loc)` avoids type-based
  // lookup entirely, preventing the exactly_once static assertion failure.
  //
  // 2. In List Mode (`UseVector` = false): Index 0 unambiguously maps to Vertex_descriptor (iterator)
  // and Index 1 maps to Edge_descriptor (iterator), making index-based access completely generic and safe across
  // both storage modes.

  // Point already exists; return existing vertex.
  if (loc.index() == 0) return std::get<0>(loc);

  auto e = std::get<1>(loc);

  // A point on a closed boundary becomes the boundary vertex: it bounds the
  // edge that contains it, and no new edge is created.  This also covers an
  // empty arrangement.
  if constexpr (Aoc::left_is_closed || Aoc::right_is_closed) {
    auto ps = arr.geometry_traits_1().parameter_space_in_x_1_object()(p);
    if constexpr (Aoc::left_is_closed) {
      if (ps == LEFT_BOUNDARY) return arr.insert_left_boundary_vertex(p);
    }
    if constexpr (Aoc::right_is_closed) {
      if (ps == RIGHT_BOUNDARY) return arr.insert_right_boundary_vertex(e, p);
    }
  }

  if (! arr.has_left_vertex(e) && ! arr.has_right_vertex(e)) return arr.insert_empty(p);

  if constexpr (Aoc::is_identified) {
    if (arr.geometry_traits_1().is_on_x_identification_1_object()(p)) return arr.insert_identification_vertex(p);
    return arr.split_edge(e, p);
  }

  // Interior edge (v_left, v_right): split it.
  return arr.split_edge(e, p);
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
