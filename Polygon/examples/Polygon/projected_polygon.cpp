#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

int main()
{
  Point_3 points[4] = { Point_3(0,1,1), Point_3(0,2,1), Point_3(0,2,2), Point_3(0,1,2) };
  bool b =  CGAL::is_simple_2(points,
                              points+4,
                              CGAL::Projection_traits_yz_3<K>());
  if (!b){
    std::cerr << "Error polygon is not simple" << std::endl;
    return 1;
  }
  return 0;
}
