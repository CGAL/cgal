#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Timer.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_3<K> DT;
typedef K::Point_3 Point_3;
typedef CGAL::Timer Timer;
typedef CGAL::Creator_uniform_3<double,Point_3>  Creator;

int main(int,char** argv)
{
  int n=atoi(argv[1]);
  std::vector<Point_3> points;
  points.reserve( n );

  CGAL::Random rng(0);
  CGAL::Random_points_in_sphere_3<Point_3,Creator> g( 1,rng);
  CGAL::cpp11::copy_n( g, n, std::back_inserter(points));

  Timer timer;
  timer.start();
  DT dt;
  dt.insert(points.begin(), points.end());
  timer.stop();
  std::size_t N = dt.number_of_vertices();
  
  std::cerr << N << std::endl;
  std::cout << timer.time() << " seconds (last"
#ifdef CGAL_TDS_USE_RECURSIVE_CREATE_STAR_3
  "CGAL_TDS_USE_RECURSIVE_CREATE_STAR_3"
#endif
#ifdef ONLY_NON_RECURSIVE
  "ONLY_NON_RECURSIVE"
#endif
#ifdef RECURSIVE_FOLLOWED_BY_UNRECURSIVE
  "RECURSIVE_FOLLOWED_BY_UNRECURSIVE"
#endif
#ifdef UNRECURSIVE_WITH_LOCAL_STACK
  "UNRECURSIVE_WITH_LOCAL_STACK"
#endif
#ifdef RECURSIVE_FOLLOWED_BY_UNRECURSIVE_WITH_LOCAL_STACK
  "RECURSIVE_FOLLOWED_BY_UNRECURSIVE_WITH_LOCAL_STACK"
#endif
  << ")\n";
  return 0;
}
