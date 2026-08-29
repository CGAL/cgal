//! \file arr_on_algebraic_segments.cpp
// Constructing a 1D arrangement on an algebraic curve using
// CGAL::Arr_algebraic_segment_traits_2 as the underlying 2D geometry kernel.

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <variant>

#include <CGAL/config.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/Geom_traits_2_adaptor_1.h>

// Setup standard 2D algebraic traits using CORE tracking
using Integer = CORE::BigInt;
using Rational = CORE::BigRat;

using Geom_traits_2 = CGAL::Arr_algebraic_segment_traits_2<Integer>;
using Polynomial_2 = Geom_traits_2::Polynomial_2;
using Algebraic_real = Geom_traits_2::Algebraic_real_1;
using Curve_2 = Geom_traits_2::Curve_2;
using X_monotone_curve_2 = Geom_traits_2::X_monotone_curve_2;

// Wrap into our 1D arrangement architecture
using Geom_traits_1 = CGAL::Arrangement_on_curve_1::Geom_traits_2_adaptor_1<Geom_traits_2>;
using Point_1 = Geom_traits_1::Point_1; // Maps to Geom_traits_2::Point_2
using Topo_traits = CGAL::Arrangement_on_curve_1::Unbounded_topology_traits<Point_1, std::string, void>;
using Arrangement = CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<Geom_traits_1, Topo_traits>;

int main() {
  // Initialize the 2D shared traits pointer
  auto traits_2_ptr = std::make_shared<const Geom_traits_2>();

  auto ctr_cv = traits_2_ptr->construct_curve_2_object();
  auto ctr_pt = traits_2_ptr->construct_point_2_object();
  auto make_xmon = traits_2_ptr->make_x_monotone_2_object();

  // --------------------------------------------------------------------------
  // 1. Construct a closed algebraic curve: the circle x^2 + y^2 - 25 = 0
  //    (center the origin, radius 5).
  // --------------------------------------------------------------------------
  Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);

  Curve_2 circle = ctr_cv(CGAL::ipower(x, 2) + CGAL::ipower(y, 2) - 25);

  // Decompose the circle into its x-monotone arcs. For a circle this yields
  // exactly two arcs: the lower semicircle (arc number 0, y <= 0) and the
  // upper semicircle (arc number 1, y >= 0), when shooting a vertical ray
  // upward from -infinity at any interior x-coordinate.
  using Make_x_monotone_result = std::variant<Point_1, X_monotone_curve_2>;
  std::vector<Make_x_monotone_result> pre_segs;
  make_xmon(circle, std::back_inserter(pre_segs));

  std::vector<X_monotone_curve_2> arcs;
  for (const auto& obj : pre_segs) {
    auto* arc_p = std::get_if<X_monotone_curve_2>(&obj);
    if (arc_p) arcs.push_back(*arc_p);
  }

  if (arcs.empty()) {
    std::cerr << "Failed to decompose the circle into x-monotone arcs.\n";
    return 1;
  }

  std::cout << "Constructed the circle x^2 + y^2 - 25 = 0 (radius 5,\n"
            << "centered at the origin), decomposed into " << arcs.size()
            << " x-monotone arc(s).\n\n";

  // --------------------------------------------------------------------------
  // 2. Wrap one of the x-monotone arcs as the "master arc" of the 1D
  //    geometry traits, and create the (initially empty) 1D arrangement.
  //
  //    Note: Compare_x_1 (defined in Geom_traits_2_adaptor_1) orders points
  //    globally via Geom_traits_2::compare_x_2_object() /
  //    compare_xy_2_object(), independent of which arc the points were
  //    constructed on. Hence any x-monotone arc of the curve is a valid
  //    master arc for the purposes of ordering points along the curve.
  // --------------------------------------------------------------------------
  Geom_traits_1 geom_traits_1(traits_2_ptr);
  Arrangement arr(std::make_shared<const Geom_traits_1>(geom_traits_1));

  // --------------------------------------------------------------------------
  // 3. Construct several points on the upper semicircle (arc number 1,
  //    i.e. y >= 0, satisfying x^2 + y^2 = 25) and insert them into the
  //    1D arrangement. The points are listed in an arbitrary order;
  //    insertion sorts them by x-coordinate along the curve.
  // --------------------------------------------------------------------------
  struct Named_point { int x; int arcno; const char* label; };
  std::vector<Named_point> samples = {
    { 3,  1, "P(x=3)"  },
    { -5, 0, "P(x=-5)" }, // left endpoint: vertical tangent, only arc 0 exists
    { 0,  1, "P(x=0)"  }, // apex of the upper arc
    { -3, 1, "P(x=-3)" },
    { 5,  0, "P(x=5)"  }, // right endpoint: vertical tangent, only arc 0 exists
    { 4,  1, "P(x=4)"  },
    { -4, 1, "P(x=-4)" },
  };

  std::cout << "Inserting points on the upper semicircle (in this order):\n";
  for (const auto& s : samples)
    std::cout << "  " << s.label << "\n";
  std::cout << "\n";

  auto v_data_map = arr.vertex_data_map();
  for (const auto& s : samples) {
    Point_1 p = ctr_pt(Algebraic_real(s.x), circle, s.arcno);
    auto v = CGAL::Arrangement_on_curve_1::insert(arr, p);
    put(v_data_map, v, std::string(s.label));
  }

  // --------------------------------------------------------------------------
  // 4. Traverse the arrangement and print the resulting sorted sequence,
  //    together with the bounded/unbounded structure of the edges.
  //
  //    Note: vertices are NOT stored in a randomly-ordered container, but the
  //    public vertices()/edges() ranges make no ordering guarantee either.
  //    To traverse the arrangement in sorted (left-to-right) order, we walk
  //    the topological linked structure explicitly, starting from the
  //    leftmost (unbounded) edge.
  // --------------------------------------------------------------------------
  std::cout << "Resulting arrangement:\n";
  std::cout << "  Number of vertices: " << arr.number_of_vertices() << "\n";
  std::cout << "  Number of edges:    " << arr.number_of_edges()    << "\n\n";

  std::cout << "Sorted vertex sequence (left to right by x-coordinate):\n";
  auto e = arr.unbounded_left_edge();
  while (arr.has_right_vertex(e)) {
    auto v = arr.right_vertex(e);
    std::cout << "  " << get(v_data_map, v) << "\n";
    e = arr.right_edge(v);
  }

  std::cout << "\nEdge structure (left to right):\n";
  size_t idx = 0;
  e = arr.unbounded_left_edge();
  while (true) {
    std::cout << "  Edge " << idx++ << ": ";
    if (arr.has_left_vertex(e)) std::cout << "[" << get(v_data_map, arr.left_vertex(e)) << "]";
    else std::cout << "(curve start)";
    std::cout << " -- ";

    if (arr.has_right_vertex(e)) std::cout << "[" << get(v_data_map, arr.right_vertex(e)) << "]";
    else std::cout << "(curve end)";
    std::cout << "\n";

    if (! arr.has_right_vertex(e)) break;
    e = arr.right_edge(arr.right_vertex(e));
  }

  // --------------------------------------------------------------------------
  // 5. Locate a query point and report whether it coincides with a vertex
  //    or falls in the interior of an edge.
  // --------------------------------------------------------------------------
  Point_1 query = ctr_pt(Algebraic_real(0), circle, 1); // the apex, previously inserted
  auto loc = CGAL::Arrangement_on_curve_1::locate(arr, query);

  std::cout << "\nLocating the apex point (x = 0):\n";
  if (loc.index() == 0) {
    auto v = std::get<0>(loc);
    std::cout << "  Found existing vertex: " << get(v_data_map, v) << "\n";
  }
  else {
    std::cout << "  Falls inside an edge (not on an existing vertex).\n";
  }

  return 0;
}
