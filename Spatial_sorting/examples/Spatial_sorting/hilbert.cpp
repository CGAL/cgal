#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/hilbert_sort.h>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef CGAL::Creator_uniform_2<double,Point>               Creator;

int main ()
{
  std::size_t size = 16;
  std::vector<Point> v; v.reserve(size);
  CGAL::points_on_square_grid_2(3.0, size,            // generate points
                                std::back_inserter(v), Creator());

  CGAL::hilbert_sort (v.begin(), v.end());            // sort

  for(std::size_t i=0; i<size; ++i)std::cout<<v[i]<<std::endl;//output
  return 0;
}
