#include <iostream>
#include <fstream>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair_2 = CGAL::Polygon_repair_2::Polygon_repair_2<Kernel>;

int main(int argc, char* argv[]) {
  // std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);

  // Square
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0))");
  // Polygon_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Bowtie
  // std::istringstream iss("POLYGON((0 0,1 1,1 0,0 1,0 0))");
  // Polygon_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Overlapping edge
  // std::istringstream iss("MULTIPOLYGON(((0 0,1 0,1 1,0 1,0 0)),((1 0,2 0,2 1,1 1,1 0)))");
  // Multipolygon_with_holes_2 p;
  // CGAL::IO::read_multi_polygon_WKT(iss, p);

  // Edge partly overlapping (start)
  // std::istringstream iss("MULTIPOLYGON(((0 0,1 0,1 1,0 1,0 0)),((1 0,2 0,2 0.5,1 0.5,1 0)))");
  // Multipolygon_with_holes_2 p;
  // CGAL::IO::read_multi_polygon_WKT(iss, p);

  // Edge partly overlapping (middle)
  // std::istringstream iss("MULTIPOLYGON(((0 0,1 0,1 1,0 1,0 0)),((1 0.25,2 0.25,2 0.75,1 0.75,1 0)))");
  // Multipolygon_with_holes_2 p;
  // CGAL::IO::read_multi_polygon_WKT(iss, p);

  // Square with hole
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0),(0.25 0.25,0.75 0.25,0.75 0.75,0.25 0.75,0.25 0.25))");
  // Polygon_with_holes_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole touching boundary
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0),(0.25 0.25,0.75 0.25,0.75 0.75,0 1,0.25 0.25))");
  // Polygon_with_holes_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole touching boundary twice
  std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0),(0.25 0.25,1 0,0.75 0.75,0 1,0.25 0.25))");
  Polygon_with_holes_2 p;
  CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole touching boundary (self-intersecting loop)
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0.25 0.25,0.75 0.25,0.75 0.75,0 1,0 0))");
  // Polygon_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole using bridge edge
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0.25 0.75,0.75 0.75,0.75 0.25,0.25 0.25,0.25 0.75,0 1,0 0))");
  // Polygon_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole outside
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0),(1.25 1.25,1.75 1.25,1.75 1.75,1.25 1.75,1.25 1.25))");
  // Polygon_with_holes_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole overlapping outer boundary
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0),(0.5 0.5,1 0.5,1 1,0.5 1,0.5 0.5))");
  // Polygon_with_holes_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  // Square with hole partly outside
  // std::istringstream iss("POLYGON((0 0,1 0,1 1,0 1,0 0),(0.75 0.75,1.25 0.75,1.25 1.25,0.75 1.25,0.75 0.75))");
  // Polygon_with_holes_2 p;
  // CGAL::IO::read_polygon_WKT(iss, p);

  CGAL::draw(p);
  Multipolygon_with_holes_2 rmp = CGAL::Polygon_repair_2::repair(p);
  CGAL::IO::write_multi_polygon_WKT(std::cout, rmp);
  CGAL::draw(rmp);

  // std::cout << "Orientations good?" << std::endl;
  // for (auto const& polygon: rmp.polygons()) {
  //   std::cout << (polygon.outer_boundary().orientation() == CGAL::COUNTERCLOCKWISE) << std::endl;
  //   for (auto const &hole: polygon.holes()) {
  //     std::cout << (hole.orientation() == CGAL::CLOCKWISE) << std::endl;
  //   }
  // }
  

  return 0;
}
