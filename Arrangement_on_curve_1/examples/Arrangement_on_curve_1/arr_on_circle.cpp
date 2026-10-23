/*! \file arr_on_circle.cpp
 * This example demonstrates `Arrangement_on_curve_1` with `Circle_traits_1`,
 * the kernel-parameterized geometry traits for points on \f$\mathbb{S}^1\f$.
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
  // The default identification point is the positive x-axis direction (1, 0).
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
    Direction_2(1, 0),    //   0 degrees (right edge normal)
    Direction_2(0, 1),    //  90 degrees (top edge normal)
    Direction_2(-1, 0),   // 180 degrees (left edge normal)
    Direction_2(0, -1),   // 270 degrees (bottom edge normal)
  };

  Arr arr_a(traits);
  for (const auto& d : normals_A) Aoc::insert(arr_a, d);

  print_circle_arr(arr_a, "square normals (A)");

  // -----------------------------------------------------------------------
  // Build arrangement B: outward edge-normal directions of a regular hexagon.
  // A hexagon has 6 edges; its normals are at 0, 60, 120, 180, 240, and 300
  // degrees.  The exact normal at 60 degrees is (1, sqrt(3)); for exact hexagon
  // geometry use a kernel that supports square roots.  Here we use integer
  // approximations for illustration:
  //     0 degrees -> (2, 0), the identification point
  //    60 degrees -> (1, 2), at about 63 degrees
  //   120 degrees -> (-1, 2), at about 117 degrees
  //   180 degrees -> (-2, 0)
  //   240 degrees -> (-1, -2), at about 243 degrees
  //   300 degrees -> (1, -2), at about 297 degrees
  //
  // Note: `Circle_traits_1` works correctly for any non-zero integer
  // `Direction_2`---no need for unit-length vectors.
  // -----------------------------------------------------------------------
  std::vector<Direction_2> normals_B = {
    Direction_2(2, 0),     //   0 degrees
    Direction_2(1, 2),     //  63 degrees
    Direction_2(-1, 2),    // 117 degrees
    Direction_2(-2, 0),    // 180 degrees
    Direction_2(-1, -2),   // 243 degrees
    Direction_2(1, -2),    // 297 degrees
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
  // Locate a query direction: (1, 1), at 45 degrees, lies between (2, 0) at 0
  // degrees and (1, 2) at about 63 degrees.
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
  // Demonstrate Is_on_x_identification_1.
  // -----------------------------------------------------------------------
  auto is_on_id = traits->is_on_x_identification_1_object();
  Direction_2 id_dir(1, 0);
  Direction_2 other_dir(0, 1);
  std::cout << "\nIs_on_x_identification_1 checks:\n";
  std::cout << "  (1,0) -> " << (is_on_id(id_dir) ? "true (identification point)" : "false") << "\n";
  std::cout << "  (0,1) -> " << (is_on_id(other_dir) ? "true (identification point)" : "false") << "\n";
  std::cout << "  (2,0) -> " << (is_on_id(Direction_2(2, 0)) ? "true (identification point)" : "false")
            << "  (same direction as (1,0), so also the identification point)\n";

  return 0;
}
