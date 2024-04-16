#include <iostream>
#include <fstream>
#include <deque>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;

int main(int argc, char* argv[])
{
  std::ifstream in((argc > 1) ? argv[1] : "data/flat.wkt");

  typedef std::vector<Point_2> MultiPoint;
  typedef std::vector<Point_2> LineString;
  typedef std::deque<LineString> MultiLineString;

  MultiPoint points;
  MultiLineString polylines;
  Multipolygon_with_holes_2 polygons;
  CGAL::IO::read_WKT(in, points,polylines,polygons);

  Multipolygon_with_holes_2 mp = CGAL::Polygon_repair::repair(polygons);
  CGAL::IO::write_multi_polygon_WKT(std::cout, mp);

  return 0;
}
