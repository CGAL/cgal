#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main(int argc, char* argv[]) {

  const char* filename = (argc > 1) ? argv[1] : "data/bridge-edge.wkt";

  std::ifstream in(filename);
  Polygon_with_holes_2 pin;
  CGAL::IO::read_polygon_WKT(in, pin);

  Multipolygon_with_holes_2 mp = CGAL::Polygon_repair::repair(pin);
  if (mp.number_of_polygons_with_holes() > 1) {
    CGAL::IO::write_multi_polygon_WKT(std::cout, mp);
  } else {
    CGAL::IO::write_polygon_WKT(std::cout, mp.polygons_with_holes()[0]);
  }

  return 0;
}
