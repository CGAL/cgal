//! \file arr_on_circle_segment.cpp
// Constructing a 1D arrangement on the upper semicircle x^2 + y^2 = 25 using
// CGAL::Arr_circle_segment_traits_2 as the underlying 2D geometry traits, and
// drawing the corresponding 2D arrangement of the half circle.

#include <iostream>
#include <string>
#include <memory>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/draw_arrangement_2.h>

#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Closed_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
#include <CGAL/Arrangement_on_curve_1/locate.h>
#include <CGAL/Arrangement_on_curve_1/Closed_geom_traits_2_adaptor_1.h>

// Kernel: Arr_circle_segment_traits_2 requires the kernel's FT to be an
// exact rational type so that circle centers and squared radii are exact.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = Kernel::FT;

// --------------------------------------------------------------------------
// 2D traits for circular arcs and line segments.
// --------------------------------------------------------------------------
using Geom_traits_2 = CGAL::Arr_circle_segment_traits_2<Kernel>;
using CoordNT = Geom_traits_2::CoordNT;
using Point_2 = Geom_traits_2::Point_2;
using X_monotone_curve_2 = Geom_traits_2::X_monotone_curve_2;
using Rat_point_2 = Geom_traits_2::Rational_point_2;
using Circle_2 = Geom_traits_2::Rational_circle_2;
using Arrangement_2 = CGAL::Arrangement_2<Geom_traits_2>;

// 1D arrangement on the upper semicircular arc.
namespace Aoc = CGAL::Arrangement_on_curve_1;

using Geom_traits_1 = Aoc::Closed_geom_traits_2_adaptor_1<Geom_traits_2>;
using Point_1 = Geom_traits_1::Point_1;   // same as Point_2 above
using Topo_traits = Aoc::Closed_topology_traits<Geom_traits_1, std::string, void>;
using Arrangement_1 = Aoc::Arrangement_on_curve_1<Geom_traits_1, Topo_traits>;

/*! Helper: print a Point_2 using double approximations of its CoordNT coords.
 */
static void print_point(const Point_2& p)
{ std::cout << "(" << CGAL::to_double(p.x()) << ", " << CGAL::to_double(p.y()) << ")"; }

/* Helper: summarise a 2D arrangement.
 */
static void print_arrangement_2(const Arrangement_2& arr) {
  std::cout << "  Vertices: " << arr.number_of_vertices()
            << "   Edges: "   << arr.number_of_edges()
            << "   Faces: "   << arr.number_of_faces()  << "\n";
  std::cout << "  Vertex coordinates (approximate):\n";
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    std::cout << "    ";
    print_point(vit->point());
    std::cout << "\n";
  }
}

int main() {
  /* Define the full circle x^2 + y^2 = 25 (center origin, radius 5, squared radius 25).
   * The two points where the circle meets the x-axis have rational coordinates: L = (-5, 0) and R = (5, 0).
   */
  Rat_point_2 center(FT(0), FT(0));
  FT sqr_radius(25);
  Circle_2 circle(center, sqr_radius, CGAL::COUNTERCLOCKWISE);

  // The endpoints (-5, 0) and (5, 0) have rational coordinates, so they
  // can be constructed with the rational Point_2 constructor.
  Point_2 pt_left (FT(-5), FT(0));   // left  x-axis crossing: (-5, 0)
  Point_2 pt_right(FT( 5), FT(0));   // right x-axis crossing: ( 5, 0)

  // Upper semicircle: counterclockwise arc from (-5,0) to (5,0) along y >= 0.
  X_monotone_curve_2 upper_arc(circle, pt_right, pt_left, CGAL::COUNTERCLOCKWISE);

  std::cout << "Circle x^2 + y^2 = 25 (center origin, radius 5).\n"
            << "Master arc: upper semicircle from (-5,0) to (5,0), y >= 0.\n\n";

  /* Build the 2D arrangement of the upper semicircle.  It has 2 vertices (the
   * x-axis crossings), 1 edge, and 1 face; the points inserted into the 1D
   * arrangement below are also inserted into it.
   */
  auto traits_2_ptr = std::make_shared<const Geom_traits_2>();
  Arrangement_2 arr_2d(traits_2_ptr.get());
  CGAL::insert(arr_2d, upper_arc);

  // Build the 1D arrangement on the upper semicircular arc.
  Geom_traits_1 geom_traits_1(traits_2_ptr, pt_left, pt_right);
  Arrangement_1 arr_1d(std::make_shared<const Geom_traits_1>(geom_traits_1));

  /* Construct several points on the upper semicircle and insert them into the 1D arrangement.
   *
   * All interior sample points satisfy x^2 + y^2 = 25, y > 0. Their x-coordinates are integers,
   * so y = sqrt(25 - x^2).
   * CoordNT = Sqrt_extension<FT, FT, ...> represents a0 + a1*sqrt(r):
   *   x-coord: pure rational ->  CoordNT(FT(ix), FT(0), FT(1))
   *   y-coord: 0 + 1*sqrt(25-ix^2) -> CoordNT(FT(0), FT(1), FT(25 - ix*ix))
   * The two endpoints where y = 0 have rational coordinates, constructed directly with Point_2(FT, FT).
   */
  auto make_upper_point = [](int ix) -> Point_2 {
    FT root(25 - ix * ix);
    CoordNT cx(FT(ix), FT(0), FT(1));   //  ix + 0*sqrt(...)  = ix
    CoordNT cy(FT(0),  FT(1), root);    //  0  + 1*sqrt(root) = sqrt(25-ix^2)
    return Point_2(cx, cy);
  };

  struct Named_point { int x; bool endpoint; const char* label; };
  std::vector<Named_point> samples = {
    {  3, false, "P(x=3)"  },
    {  0, false, "P(x=0)"  }, // apex  of the upper arc, y = 5
    { -5, true,  "P(x=-5)" }, // left  endpoint of the upper arc, y = 0
    { -3, false, "P(x=-3)" },
    {  5, true,  "P(x=5)"  }, // right endpoint of the upper arc, y = 0
    {  4, false, "P(x=4)"  },
    { -4, false, "P(x=-4)" },
  };

  std::cout << "Inserting points on the upper semicircle (in this order):\n";
  for (const auto& s : samples) std::cout << "  " << s.label << "\n";
  std::cout << "\n";

  auto v_data_map = arr_1d.vertex_data_map();
  for (const auto& s : samples) {
    // exact rational endpoint or algebraic interior point
    Point_1 p = s.endpoint ? Point_2(FT(s.x), FT(0)) : make_upper_point(s.x);
    auto v = CGAL::Arrangement_on_curve_1::insert(arr_1d, p);
    CGAL::insert_point(arr_2d, p);
    put(v_data_map, v, std::string(s.label));
  }

  std::cout << "Resulting 1D arrangement:\n";
  std::cout << "  Number of vertices: " << arr_1d.number_of_vertices() << "\n";
  std::cout << "  Number of edges:    " << arr_1d.number_of_edges()    << "\n\n";

  /* Traverse the 1D arrangement in sorted (left-to-right) order, starting at
   * the leftmost edge.  The arrangement has number_of_edges() edges, and every
   * edge but the last has a right vertex.
   */
  auto v_pnt_map = arr_1d.vertex_point_map();
  auto print_vertex = [&](Arrangement_1::Vertex_const_descriptor v) {
    std::cout << "  " << get(v_data_map, v) << "  ~  ";
    print_point(get(v_pnt_map, v));
    std::cout << "\n";
  };
  auto next = [&](auto e) { return arr_1d.right_edge(arr_1d.right_vertex(e)); };

  std::cout << "Sorted vertex sequence (left to right by x-coordinate):\n";
  auto e = arr_1d.left_boundary_edge();
  if (arr_1d.has_left_vertex(e)) print_vertex(arr_1d.left_vertex(e));
  for (std::size_t i = 0; i < arr_1d.number_of_edges(); ++i) {
    if (i > 0) e = next(e);
    if (arr_1d.has_right_vertex(e)) print_vertex(arr_1d.right_vertex(e));
  }

  std::cout << "\nEdge structure (left to right):\n";
  e = arr_1d.left_boundary_edge();
  for (std::size_t i = 0; i < arr_1d.number_of_edges(); ++i) {
    if (i > 0) e = next(e);
    std::cout << "  Edge " << i << ": ";
    if (arr_1d.has_left_vertex(e)) std::cout << "[" << get(v_data_map, arr_1d.left_vertex(e)) << "]";
    else std::cout << "(arc start)";
    std::cout << " -- ";
    if (arr_1d.has_right_vertex(e)) std::cout << "[" << get(v_data_map, arr_1d.right_vertex(e)) << "]";
    else std::cout << "(arc end)";
    std::cout << "\n";
  }

  /* Locate a query point and report whether it coincides with a vertex.
   */
  Point_1 query = make_upper_point(0);   // the apex (0, 5), already inserted
  auto loc = CGAL::Arrangement_on_curve_1::locate(arr_1d, query);

  std::cout << "\nLocating the apex point (x = 0, y = 5):\n";
  if (loc.index() == 0) {
    auto v = std::get<0>(loc);
    std::cout << "  Found existing vertex: " << get(v_data_map, v) << "\n";
  }
  else {
    std::cout << "  Falls inside an edge (not on an existing vertex).\n";
  }

  std::cout << "2D arrangement of the half circle:\n";
  print_arrangement_2(arr_2d);
  std::cout << "\n";

  CGAL::Graphics_scene_options<Arrangement_2, typename Arrangement_2::Vertex_const_handle,
                               typename Arrangement_2::Halfedge_const_handle,
                               typename Arrangement_2::Face_const_handle> gso;
  gso.ignore_all_vertices(false);
  gso.colored_vertex = [](const Arrangement_2&, Arrangement_2::Vertex_const_handle) -> bool { return true; };
  gso.vertex_color = [](const Arrangement_2&, Arrangement_2::Vertex_const_handle) -> CGAL::IO::Color
  { return CGAL::IO::Color(128, 0, 128); };
  gso.colored_face = [](const Arrangement_2&, Arrangement_2::Face_const_handle) -> bool { return false; };
  CGAL::draw(arr_2d, gso, "circle");

  return 0;
}
