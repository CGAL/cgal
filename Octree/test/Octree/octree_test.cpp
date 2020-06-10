
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Kernel, Point_set, typename Point_set::Point_map, typename Point_set::Vector_map>
        Octree;

int main(void) {

  Point_set points;
  points.insert({0, 0, 1});

  auto point_map = points.point_map();
  auto normal_map = points.normal_map();

  Octree octree(points, point_map, normal_map);

  std::cout << *octree.root();

  return 0;
}