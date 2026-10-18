// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

/*! \file arr_on_circle.cpp
 * Example: 1D arrangement on a circle (identified parameter space).
 *
 * This example demonstrates `Arrangement_on_curve_1` with `Circle_traits_1`,
 * the kernel-parameterised geometry traits for points on \f$\mathbb{S}_1\f$.
 * Points are represented as `Direction_2` objects from the kernel.
 */

#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_curve_1/Circle_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Identified_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/locate.h>
#include <CGAL/Arrangement_on_curve_1/overlay.h>
#include <CGAL/Arrangement_on_curve_1/Default_overlay_visitor.h>

// ---------------------------------------------------------------------------
// Type setup
// ---------------------------------------------------------------------------
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits_1 = CGAL::Arrangement_on_curve_1::Circle_traits_1<Kernel>;
using Direction_2 = Kernel::Direction_2;
using Point_1 = Traits_1::Point_1;
using Topol_traits = CGAL::Arrangement_on_curve_1::Identified_topology_traits<Traits_1>;
using Arr = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Traits_1, Topol_traits>;

namespace Aoc = CGAL::Arrangement_on_curve_1;

// ---------------------------------------------------------------------------
// Helper: pretty-print a Direction_2 as "dx/dy"
// ---------------------------------------------------------------------------
static std::ostream& print_dir(std::ostream& os, const Direction_2& d)
{ return os << "(" << d.dx() << ", " << d.dy() << ")"; }

// ---------------------------------------------------------------------------
// Pretty-print an arrangement on a circle by walking the cyclic structure.
// ---------------------------------------------------------------------------
static void print_circle_arr(const Arr& arr, const std::string& name) {
  std::cout << "Arrangement on circle: " << name << "\n";
  std::cout << "  vertices = " << arr.number_of_vertices() << ", edges = "    << arr.number_of_edges() << "\n";

  if (arr.empty()) { std::cout << "  (empty)\n\n"; return; }

  auto& topo = arr.topology_traits();
  auto  pmap = arr.vertex_point_map();

  auto boundary_e = topo.left_boundary_edge();
  auto cur_e = boundary_e;
  std::cout << "  cyclic vertex sequence (direction dx, dy):\n";
  std::size_t count = 0;
  do {
    auto v = topo.right_vertex(cur_e);
    std::cout << "    ";
    print_dir(std::cout, get(pmap, v));
    std::cout << "\n";
    cur_e = topo.right_edge(v);
    ++count;
  } while (cur_e != boundary_e && count <= arr.number_of_vertices());
  std::cout << "\n";
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
  // The default seam is the positive x-axis direction (1, 0).
  auto traits = std::make_shared<const Traits_1>();

  // -----------------------------------------------------------------------
  // Build arrangement A: outward edge-normal directions of a square.
  // A square has 4 edges; normals are the four axis-aligned unit directions.
  //   right (1,0) the identification point.
  //   top (0,1)
  //   left (-1,0)
  //   bottom (0,-1)
  // -----------------------------------------------------------------------
  std::vector<Direction_2> normals_A = {
    Direction_2(1,0),     // \f$0^\circ\f$  (right edge normal)
    Direction_2(0, 1),    // \f$90^\circ\f$  (top edge normal)
    Direction_2(-1, 0),   // \f$180^\circ\f$ (left edge normal)
    Direction_2(0, -1),   // \f$270^\circ\f$ (bottom edge normal)
  };

  Arr arr_a(traits);
  for (const auto& d : normals_A) Aoc::insert(arr_a, d);

  print_circle_arr(arr_a, "square normals (A)");

  // -----------------------------------------------------------------------
  // Build arrangement B: outward edge-normal directions of a regular hexagon.
  // A hexagon has 6 edges; normals at \f$0^\circ\f$, \f$60^\circ\f$,
  // \f$120^\circ\f$, \f$180^\circ\f$, \f$240^\circ\f$, \f$300^\circ\f$.  Using
  // exact integer `Direction_2` vectors:
  //   \f$0^\circ  \rightarrow (2, 0)\f$-- the identification point
  //   \f$60^\circ \rightarrow (1, \sqrt{3})\f$---approximate: (1, 2) is close but not exact;
  //                     for exact hex geometry use a sqrt extension kernel.
  //                     Here we use integer approximations for illustration.
  //  \f$120^\circ \rightarrow (-1, 2)\f$
  //  \f$180^\circ \rightarrow (-2, 0)\f$
  //  \f$240^\circ \rightarrow (-1, -2)\f$
  //  \f$300^\circ \rightarrow (1, -2)\f$
  //
  // Note: `Circle_traits_1` works correctly for any non-zero integer
  // `Direction_2`---no need for unit-length vectors.
  // -----------------------------------------------------------------------
  std::vector<Direction_2> normals_B = {
    Direction_2(2, 0),     // \f$\approx   0^\circ\f$
    Direction_2(1, 2),     // \f$\approx  60^\circ\f$
    Direction_2(-1, 2),    // \f$\approx 120^\circ\f$
    Direction_2(-2, 0),    // \f$\approx 180^\circ\f$
    Direction_2(-1, -2),   // \f$\approx 240^\circ\f$
    Direction_2( 1, -2),   // \f$\approx 300^\circ\f$
  };

  Arr arr_b(traits);
  for (const auto& d : normals_B) Aoc::insert(arr_b, d);

  print_circle_arr(arr_b, "hexagon normals (B)");
  // -----------------------------------------------------------------------
  // Overlay: the result arrangement contains a vertex for every direction
  // that is a normal in A, in B, or in both.
  // -----------------------------------------------------------------------
  Arr arr_res(traits);
  using DefaultVisitor = Aoc::Default_overlay_visitor<Arr, Arr, Arr>;
  DefaultVisitor visitor(arr_a, arr_b, arr_res);
  Aoc::overlay(arr_a, arr_b, arr_res, visitor);

  print_circle_arr(arr_res, "overlay A and B");

  // -----------------------------------------------------------------------
  // Locate a query direction: (1, 1) ≈ 45°, between 0° and ≈60°.
  // -----------------------------------------------------------------------
  Direction_2 query(1, 1);
  auto loc = Aoc::locate(arr_res, query);
  std::cout << "locate";
  print_dir(std::cout, query);
  std::cout << " in arr_res:\n";
  if (loc.index() == 0) {
    auto v = std::get<0>(loc);
    std::cout << "  -> coincides with vertex at ";
    print_dir(std::cout, get(arr_res.vertex_point_map(), v));
    std::cout << "\n";
  }
  else std::cout << "  -> falls inside an edge\n";

  // -----------------------------------------------------------------------
  // Demonstrate Is_on_identification_1.
  // -----------------------------------------------------------------------
  auto is_seam = traits->is_on_identification_1_object();
  Direction_2 seam_dir(1, 0);
  Direction_2 other_dir(0, 1);
  std::cout << "\nIs_on_identification_1 checks:\n";
  std::cout << "  (1,0) -> " << (is_seam(seam_dir)  ? "true (seam)"  : "false") << "\n";
  std::cout << "  (0,1) -> " << (is_seam(other_dir) ? "true (seam)"  : "false") << "\n";
  std::cout << "  (2,0) -> " << (is_seam(Direction_2(2,0)) ? "true (seam)" : "false")
            << "  (same direction as (1,0), so also the seam)\n";

  return 0;
}
