#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/spatial_sort.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef K::Point_3                                                Point_3;
typedef CGAL::Creator_uniform_3<double,Point_3>                   Creator_3;

int main()
{
  std::size_t pt_nb = 10000;
  std::vector<Point_3> points;
  points.reserve(pt_nb);

  CGAL::Random_points_in_cube_3<Point_3> gen(1.0);

  for(std::size_t i=0; i<pt_nb; ++i)
    points.push_back(*gen++);

  // By default sequential
  CGAL::spatial_sort(points.begin(),points.end());

  // Add the template argument to switch on parallelism if available
  // You can also use Parallel_tag if you know that TBB is enabled
  CGAL::spatial_sort<CGAL::Parallel_if_available_tag>(points.begin(), points.end());

  return 0;
}
