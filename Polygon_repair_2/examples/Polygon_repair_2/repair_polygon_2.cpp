#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_repair_2/Polygon_repair_2.h>

// #include <CGAL/IO/WKT.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Multipolygon_with_holes_2<Kernel> Multipolygon;
typedef CGAL::Polygon_repair_2::Polygon_repair_2<Kernel> Polygon_repair;

int main(int argc, char* argv[]) {
  // std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);

  Polygon_repair pr;

  return 0;
}
