#define CGAL_NO_CDT_2_WARNING

#include <iostream>
#include <fstream>
#include <sstream>

// work around for old compilers (Apple clang < 11 for example)
#define HAS_FILESYSTEM 1
#if defined(__has_include)
#if !__has_include(<filesystem>)
#undef HAS_FILESYSTEM
#define HAS_FILESYSTEM 0
#endif
#endif


#if HAS_FILESYSTEM

#include <filesystem>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main() {

  std::string in = "POLYGON((0.03 0.02,0.97 0.01,0.99 0.96,0.04 0.98,0.03 0.02),(0.5 0.5,1.5 0.5,0.5 0.5,1.5 0.7,0.5 0.5,1.5 0.9,0.5 0.5,1.5 1.1,0.5 0.5,1.5 1.3,0.5 0.5,1.5 1.5,0.5 0.5,1.3 1.5,0.5 0.5,1.1 1.5,0.5 0.5,0.9 1.5,0.5 0.5,0.7 1.5,0.5 0.5,0.5 1.5,0.5 0.5,0.3 1.5,0.5 0.5,0.1 1.5,0.5 0.5,-0.1 1.5,0.5 0.5,-0.3 1.5,0.5 0.5,-0.5 1.5,0.5 0.5,-0.5 1.3,0.5 0.5,-0.5 1.1,0.5 0.5,-0.5 0.9,0.5 0.5,-0.5 0.9,0.5 0.5,-0.5 0.7,0.5 0.5,-0.5 0.5,0.5 0.5,-0.5 0.3,0.5 0.5,-0.5 0.1,0.5 0.5,-0.5 -0.1,0.5 0.5,-0.5 -0.3,0.5 0.5,-0.5 -0.5,0.5 0.5,-0.3 -0.5,0.5 0.5,-0.1 -0.5,0.5 0.5,0.1 -0.5,0.5 0.5,0.3 -0.5,0.5 0.5,0.5 -0.5,0.5 0.5,0.7 -0.5,0.5 0.5,0.9 -0.5,0.5 0.5,1.1 -0.5,0.5 0.5,1.3 -0.5,0.5 0.5,1.5 -0.5,0.5 0.5,1.5 -0.3,0.5 0.5,1.5 -0.1,0.5 0.5,1.5 0.1,0.5 0.5,1.5 0.3,0.5 0.5))";
  std::istringstream iss(in);
  Multipolygon_with_holes_2 rmp;

  Polygon_with_holes_2 p;
  CGAL::IO::read_polygon_WKT(iss, p);
  CGAL::draw(p);
  Polygon_repair pr;
  for (auto const edge: p.outer_boundary().edges()) {
    pr.triangulation().even_odd_insert_constraint(edge.source(), edge.target());
  } int spikes = 20;
  for (auto const& hole: p.holes()) {
    for (auto const edge: hole.edges()) {
      if (spikes-- <= 0) break;
      pr.triangulation().even_odd_insert_constraint(edge.source(), edge.target());
    }
  }
  pr.label_triangulation_even_odd();
  pr.reconstruct_multipolygon();
  rmp = CGAL::Polygon_repair::repair(p, CGAL::Polygon_repair::Even_odd_rule());
  std::ostringstream oss;
  CGAL::IO::write_multi_polygon_WKT(oss, rmp);
  std::string out = oss.str();
  std::cout << "\tin:  " << in << std::endl;
  std::cout << "\tout: " << out;
  CGAL::draw(rmp);

  return 0;
}

#else

int main()
{
  std::cout << "Warning: filesystem feature is not present on the system, nothing will be tested\n";
  return 0;
}

#endif
