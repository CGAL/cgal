#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>

template <typename Kernel>
void test()
{
  typedef typename Kernel::Point_3 Point;
  typedef CGAL::Point_set_3<Point> Point_set;
  typedef CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map> Octree;

  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(100);
  for (std::size_t i = 0; i < 100; ++i)
    points.insert(*(generator++));

  Octree octree(points, points.point_map());
  octree.refine();
  octree.grade();
}

int main (int, char**)
{
  test<CGAL::Simple_cartesian<float> >();
  test<CGAL::Simple_cartesian<double> >();
  test<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test<CGAL::Exact_predicates_exact_constructions_kernel>();
  return EXIT_SUCCESS;
}
