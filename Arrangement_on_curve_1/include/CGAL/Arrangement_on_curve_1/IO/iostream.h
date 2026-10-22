// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_IO_IOSTREAM_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_IO_IOSTREAM_H

#include <cstddef>
#include <iostream>

#include <CGAL/enum.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// ============================================================================
// The format: the number of vertices followed by their points in increasing
// order along the curve, each on its own line.  In an identified arrangement,
// the identification point, if present, is last.  Points are written and read
// with the stream operators of `Point_1`.
// ============================================================================

/*! inserts the arrangement `arr` into the output stream `os`.
 */
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
std::ostream& operator<<(std::ostream& os,
                         const Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr) {
  using Aoc = Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>;

  auto pm = arr.vertex_point_map();
  os << arr.number_of_vertices() << '\n';
  if (arr.empty()) return os;

  auto e = arr.left_boundary_edge();
  if constexpr (Aoc::is_identified) {
    // Every edge has a right vertex; the walk ends at the identification
    // vertex, if it exists, since it is the left vertex of left_boundary_edge().
    for (auto n = arr.number_of_vertices(); n > 0; --n) {
      auto v = arr.right_vertex(e);
      os << get(pm, v) << '\n';
      e = arr.right_edge(v);
    }
  }
  else {
    if (arr.has_left_vertex(e)) os << get(pm, arr.left_vertex(e)) << '\n';
    while (arr.has_right_vertex(e)) {
      auto v = arr.right_vertex(e);
      os << get(pm, v) << '\n';
      // Both edges of the right boundary vertex are its left edge.
      if (arr.right_edge(v) == e) break;
      e = arr.right_edge(v);
    }
  }
  return os;
}

/*! extracts an arrangement from the input stream `is` into `arr`, replacing
 * its content.  If the input is malformed, the `failbit` of `is` is set and
 * `arr` is left empty.
 */
template <typename GeometryTraits_1, typename TopologyTraits, bool BinarySearch>
std::istream& operator>>(std::istream& is,
                         Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>& arr) {
  using Aoc = Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits, BinarySearch>;
  using Point_1 = typename Aoc::Point_1;
  using Vertex_descriptor = typename Aoc::Vertex_descriptor;

  std::size_t n;
  if (! (is >> n)) return is;
  arr.clear();

  auto fail = [&]() -> std::istream& {
    arr.clear();
    is.setstate(std::ios::failbit);
    return is;
  };

  const auto& gt = arr.geometry_traits_1();
  auto cmp = gt.compare_x_1_object();
  auto pm = arr.vertex_point_map();

  // The points are appended at the right end; e is the rightmost edge.
  auto e = arr.left_boundary_edge();
  Vertex_descriptor v_prev {};
  bool has_prev = false;

  for (std::size_t i = 0; i < n; ++i) {
    Point_1 p;
    if (! (is >> p)) return fail();
    const bool is_last = (i + 1 == n);

    if constexpr (Aoc::is_identified) {
      // Compare_x_1 must not be called on the identification point, which is
      // last; if it does not exist, left_boundary_edge() keeps the wrap-around
      // part (p, v_1) of every insertion.
      if (gt.is_on_x_identification_1_object()(p)) {
        if (! is_last) return fail();
        if (arr.empty()) arr.insert_empty(p);
        else arr.insert_identification_vertex(p);
        continue;
      }
      if (has_prev && (cmp(get(pm, v_prev), p) != SMALLER)) return fail();
      v_prev = arr.empty() ? arr.insert_empty(p) : arr.insert_after(arr.left_vertex(e), p);
    }
    else {
      if (has_prev && (cmp(get(pm, v_prev), p) != SMALLER)) return fail();
      if constexpr (Aoc::left_is_closed || Aoc::right_is_closed) {
        auto ps = gt.parameter_space_in_x_1_object()(p);
        if constexpr (Aoc::left_is_closed) {
          if (ps == Parameter_space::LEFT_BOUNDARY) {
            if (i != 0) return fail();
            v_prev = arr.insert_left_boundary_vertex(p);
            has_prev = true;
            continue;
          }
        }
        if constexpr (Aoc::right_is_closed) {
          if (ps == Parameter_space::RIGHT_BOUNDARY) {
            if (! is_last) return fail();
            v_prev = arr.insert_right_boundary_vertex(e, p);
            has_prev = true;
            continue;
          }
        }
      }
      v_prev = arr.split_edge(e, p);
      e = arr.right_edge(v_prev);
    }
    has_prev = true;
  }
  return is;
}

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
