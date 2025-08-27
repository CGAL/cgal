#include <CGAL/Polygon_repair/repair.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_repair = CGAL::Polygon_repair::Polygon_repair<Kernel>;

int main() {

  std::string polygon("POLYGON((0 0, 3 0, 3 3, 4 3,0 0))");
  std::istringstream  in(polygon);
  Polygon_with_holes_2 pwh;
  CGAL::IO::read_polygon_WKT(in, pwh);

  Polygon_2 outer = pwh.outer_boundary();

  Multipolygon_with_holes_2 rmp;

  rmp = CGAL::Polygon_repair::repair(outer, CGAL::Polygon_repair::Non_zero_rule());

  rmp = CGAL::Polygon_repair::repair(pwh, CGAL::Polygon_repair::Non_zero_rule());

  rmp = CGAL::Polygon_repair::repair(rmp, CGAL::Polygon_repair::Non_zero_rule());

std::cout << "done" << std::endl;
  return 0;
}

