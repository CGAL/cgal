#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Simple_cartesian.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;

int main(void) {

  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  points.insert({-1, 1, -1});
  points.insert({1, 1, -1});
  points.insert({-1, -1, 1});
  points.insert({1, -1, 1});
  points.insert({-1, 1, 1});
  points.insert({1, 1, 1});

  points.insert({-1, -1, -1.1});
  points.insert({-1, -1, -1.2});
  points.insert({-1, -1, -1.3});
  points.insert({-1, -1, -1.4});
  points.insert({-1, -1, -1.5});
  points.insert({-1, -1, -1.6});
  points.insert({-1, -1, -1.7});
  points.insert({-1, -1, -1.8});
  points.insert({-1, -1, -1.9});

  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  std::cout << "root: " << octree.local_coordinates(octree.root()) << std::endl;
  std::cout << "first child: " << octree.local_coordinates(octree.child(octree.root(), 0)) << std::endl;
  std::cout << "fifth child: " << octree.local_coordinates(octree.child(octree.root(), 4)) << std::endl;
  std::cout << "fifth child of first child: "
            << octree.local_coordinates(octree.child(octree.child(octree.root(), 0), 4)) << std::endl;

  return 0;
}
