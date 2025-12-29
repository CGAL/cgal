#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>
#ifdef CGAL_USE_BASIC_VIEWER
# include <CGAL/draw_polygon_2.h>
# include <CGAL/draw_polygon_with_holes_2.h>
# include <CGAL/draw_multipolygon_with_holes_2.h>
#endif

#include <iostream>
#include <sstream>
#include <string>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

// These are polygons of Figs 16.1 and 16.2 in the user manual + some customs
std::vector<std::pair<std::string, bool> > inputs =
  {{
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9)))", true},
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 5)))", true},
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)))", true},
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)), ((6 6, 8 3, 8 8)))", true},
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)), ((5 6, 8 3, 8 8)))", true},
    {"MULTIPOLYGON(((0 0, 6 0, 3 4)), ((6 0, 12 0, 9 4)), ((3 4, 9 4, 6 8)))", true}, // triforce
    {"MULTIPOLYGON(((0 0, 6 0, 3 4)), ((6 0, 12 0, 9 4)), ((3 4, 9 4, 6 8)), ((5 2, 7 2, 6 3)))", true}, // triforce + inner polygon
    {"MULTIPOLYGON(((0 0, 10 0, 10 10, 0 10), (1 9, 9 9, 9 1, 1 1)), ((2 2, 8 2, 8 8, 2 8), (3 7, 7 7, 7 3, 3 3)), ((4 4, 6 4, 6 6, 4 6)))", true}, // polygon with hole within a polygon with hole

    {"MULTIPOLYGON(((3 9, 8 9, 5 6, 7 3, 1 1)))", false}, // bad orientation
    {"MULTIPOLYGON(((1 1, 10 5, 5 4, 5 6, 10 5, 8 9, 3 9)))", false}, // touching vertices
    {"MULTIPOLYGON(((1 1, 10 5, 5 4, 5 6, 5 5, 8 9, 3 9)))", false}, // partial needle
    {"MULTIPOLYGON(((1 1, 10 5, 5 4, 5 6, 5 4, 8 9, 3 9)))", false}, // needle
    {"MULTIPOLYGON(((1 1, 7 3, 3 9, 8 9)))", false}, // edge-edge intersection
    {"MULTIPOLYGON(((1 1, 7 3, 2 5, 8 9, 3 9)))", false}, // vertex touching edge
    {"MULTIPOLYGON(((1 1, 7 3)))", false}, // combinatorially degenerate polygon
    {"MULTIPOLYGON(((1 1, 7 3, 4 2)))", false}, // geometrically degenerate polygon
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (2 2, 4 4)))", false}, // combinatorially degenerate hole
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (2 2, 4 4, 3 3)))", false}, // geometrically degenerate hole
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (5 5, 4 6, 3 4)))", false}, // bad hole orientation
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (2 2, 2 4, 4 4, 4 6)))", false}, // self-intersecting hole
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (2 2, 2 4, 3 4, 3 6, 4 6, 3 4)))", false}, // hole w/ non-manifold vertex
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (6 6, 8 3, 8 8)))", false}, // hole fully outside of hull
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (6 4, 4 6, 6 8)))", false}, // hole partially outside of hull
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (7 3, 4 6, 5 6)))", false}, // hole sharing an edge with hull
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9, 2 5), (4 4, 2 5, 5 6)))", false}, // disconnected interiors
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9, 2 5), (4 4, 3 5, 5 5), (3 4, 4 6, 5 4)))", false}, // overlapping holes
    // {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9, 2 5), (4 4, 3 5, 5 5), (3 5, 4 6, 5 5)))", false}, // holes sharing an edge
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9, 2 5), (2 2, 4 6, 6 4), (3 3, 4 5, 5 4)))", false}, // holes within a hole
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9)), ((0 0, 10 0, 10 10, 0 10)))", false}, // fully overlapping polygons
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)), ((8 3, 6 8, 8 8)))", false}, // partly overlapping polygons
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)), ((8 3, 5 6, 6 7)))", false}, // polygons partially sharing an edge
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)), ((8 3, 5 6, 8 9)))", false}, // polygons fully sharing an edge, on different sides
    {"MULTIPOLYGON(((1 1, 7 3, 5 6, 8 9, 3 9), (3 4, 4 6, 5 6)), ((5 6, 8 9, 4 8)))", false} // polygons fully sharing an edge, on the same side
  }};

bool test(const std::string& name,
          const std::string& in,
          const bool expected_result)
{
  std::cout << "\n== Testing Input " << name << " ==" << std::endl;

  std::istringstream iss(in);
  Multipolygon_with_holes_2 mph;
  if (!CGAL::IO::read_multi_polygon_WKT(iss, mph)) {
    std::cerr << "*** Error: Failed to read polygon ***" << std::endl;
    return false;
  }

#ifdef CGAL_USE_BASIC_VIEWER
  CGAL::draw(mph);
#endif

  const bool is_valid = CGAL::Polygon_repair::is_valid(mph);

  std::cout << "Input " << name << " is " << (is_valid ? "" : "NOT ") << "valid" << std::endl;

  if (is_valid != expected_result) {
    std::cerr << "*** Error: expected " << expected_result << " and got " << is_valid << " ***" << std::endl;
  }

  const bool successful_diag = (is_valid == expected_result);

  // might as well repair since we are here
  auto rmp = CGAL::Polygon_repair::repair(mph, CGAL::Polygon_repair::Even_odd_rule());

#ifdef CGAL_USE_BASIC_VIEWER
  CGAL::draw(rmp);
#endif

  const bool successful_repair = CGAL::Polygon_repair::is_valid(rmp);
  std::cout << "Repaired " << name << " is " << (successful_repair ? "" : "NOT ") << "valid" << std::endl;
  if (!successful_repair) {
    std::cerr << "*** Error: failed to repair " << (is_valid ? "valid" : "invalid") << " input ***" << std::endl;
  }

  return (successful_diag && successful_repair);
}

int main()
{
  bool success = true;
  for (std::size_t i=0; i<inputs.size(); ++i) {
    success = test(std::to_string(i), inputs[i].first, inputs[i].second) && success;
  }

  assert(success);

  std::cout << "\nDone" << std::endl;
  return EXIT_SUCCESS;
}

