
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/spatial_sort.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <vector>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Creator_uniform_3<double,Point_3> Creator_3;

int main(int argc, char* argv[])
{
  std::vector<Point_3> points;
  v.reserve (2000);

  CGAL::Random_points_in_cube_3<Point_3> gen (1.0);
  
  for (int i = 0; i < 2000; ++i){
    v.push_back (*gen++);
  }
#ifdef CGAL_LINKED_WITH_TBB
  CGAL::spatial_sort<CGAL::Parallel_tag>(points.begin(),points.end());
#else
  CGAL::spatial_sort(points.begin(),points.end());
#endif

  return 0;
}
