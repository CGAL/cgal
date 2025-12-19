#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <iostream>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

int main()
{
  Polygon_2 polygon;
  for (int i = 0; i < 200000; ++i) {
      polygon.push_back(Point_2(i* 1.04663, 0));
  }
  polygon.push_back(Point_2( 3.1415, 3.1415));

  auto ar = CGAL::polygon_area_2(polygon.vertices_begin(), polygon.vertices_end(), Kernel());

  std::cout << "done" << std::endl;
  return 0;
}