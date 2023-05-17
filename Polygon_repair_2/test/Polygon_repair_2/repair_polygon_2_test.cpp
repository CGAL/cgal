#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_repair_2/repair.h>
#include <CGAL/IO/WKT.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>                       Polygon_with_holes_2;

int main(int argc, char* argv[])
{
  std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);
  return 0;
}
