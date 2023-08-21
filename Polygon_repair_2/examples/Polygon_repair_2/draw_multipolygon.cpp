#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main(int argc, char* argv[]) {
  std::ifstream in("data/nesting.wkt");
  Multipolygon_with_holes_2 mp;
  CGAL::IO::read_multi_polygon_WKT(in, mp);
  CGAL::draw(mp);

  return 0;
}
