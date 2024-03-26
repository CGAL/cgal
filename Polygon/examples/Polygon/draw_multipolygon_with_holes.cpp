#include <iostream>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/draw_multipolygon_with_holes_2.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main() {
  std::string wkt = "MULTIPOLYGON(((0 0,1 0,1 1,0 1,0 0),(0.1 0.1,0.1 0.9,0.9 0.9,0.9 0.1,0.1 0.1)),((0.2 0.2,0.8 0.2,0.8 0.8,0.2 0.8,0.2 0.2),(0.3 0.3,0.3 0.7,0.7 0.7,0.7 0.3,0.3 0.3)))";
  std::istringstream iss(wkt);
  Multipolygon_with_holes_2 mp;
  CGAL::IO::read_multi_polygon_WKT(iss, mp);
  CGAL::draw(mp);

  return 0;
}
