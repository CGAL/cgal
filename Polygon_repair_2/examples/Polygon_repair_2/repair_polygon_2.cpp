#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main(int argc, char* argv[]) {
  std::ifstream in("data/bridge-edge.wkt");
  Polygon_with_holes_2 pin;
  CGAL::IO::read_polygon_WKT(in, pin);

  Multipolygon_with_holes_2 mp = CGAL::Polygon_repair_2::repair_odd_even(pin);
  if (mp.number_of_polygons() > 1) {
    CGAL::IO::write_multi_polygon_WKT(std::cout, mp);
  } else {
    CGAL::IO::write_polygon_WKT(std::cout, mp.polygons()[0]);
  }
  
  return 0;
}
