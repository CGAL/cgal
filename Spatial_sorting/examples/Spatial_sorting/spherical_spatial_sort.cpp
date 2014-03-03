#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/spatial_sort.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_3                                              Point;
typedef CGAL::Creator_uniform_3<double,Point>                   Creator_3;

int main ()
{
  std::size_t size = 32;
  std::vector<Point> v; v.reserve(size);

  CGAL::Random random (42);
  CGAL::Random_points_on_sphere_3<Point> gen(1.0, random);
  for (int i = 0; i < size; ++i) v.push_back(*gen++);

  CGAL::spherical_spatial_sort(v.begin(), v.end());                 // sort

  for(std::size_t i=0; i<size; ++i) std::cout << v[i] << std::endl; //output
  return 0;
}
