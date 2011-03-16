
#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
#include <CGAL/point_generators_3.h>

#include "Orientation_3.h"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
#if 1
typedef CGAL::internal::Static_filters_predicates::Orientation_3_benchmark<Point_3> Predicate;
#else

   struct Predicate {
   int operator()(const Point_3& p,const Point_3& q,const Point_3& r,const Point_3& s) const
   {
     return (int)(p < q)  + (int)(r<s);
   }
 };
#endif

int main()
{
  const int N = 100000000;
  std::vector<Point_3> points;
  points.reserve(N);
  CGAL::Random_points_in_sphere_3<Point_3> g( 100.0);
  CGAL::copy_n( g, N, std::back_inserter(points));
  CGAL::Timer timer;
  timer.start();
  int res=0;
  Predicate predicate;
  for(int i = 0; i < N-4; i++){
    res += predicate(points[i], points[i+1], points[i+2], points[i+3]);
  }  
  timer.stop();
  std::cout << res << std::endl << timer.time() << " sec" << std::endl;
}
