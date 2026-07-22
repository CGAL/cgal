#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main(int argc, char* argv[])
{
  std::ifstream in((argc > 1) ? argv[1] : CGAL::data_file_path("wkt/issue.wkt"));
  Multipolygon_with_holes_2 pin;
  CGAL::IO::read_multi_polygon_WKT(in, pin);

  Multipolygon_with_holes_2 mp = CGAL::Polygon_repair::repair(pin);
  CGAL::IO::write_multi_polygon_WKT(std::cout, mp);

  return 0;
}
