#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>
Octree;
typedef Octree::Node Node;
typedef Octree::Traits Traits;

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

  std::cout << octree << std::endl;

  // Root node should have no siblings
  assert(octree.root().adjacent_node(0).is_null());
  assert(octree.root().adjacent_node(1).is_null());
  assert(octree.root().adjacent_node(2).is_null());
  assert(octree.root().adjacent_node(3).is_null());
  assert(octree.root().adjacent_node(4).is_null());
  assert(octree.root().adjacent_node(5).is_null());

  // Left Top Front node should have siblings to the Right, Down, and Back
  auto left_top_back = octree.root()[Traits::LEFT_TOP_BACK];

  assert(octree.root()[Traits::RIGHT_TOP_BACK] == left_top_back.adjacent_node(Traits::RIGHT));
  assert(octree.root()[Traits::LEFT_BOTTOM_BACK] == left_top_back.adjacent_node(Traits::DOWN));
  assert(octree.root()[Traits::LEFT_TOP_FRONT] == left_top_back.adjacent_node(Traits::FRONT));
  assert(left_top_back.adjacent_node(Traits::LEFT).is_null());
  assert(left_top_back.adjacent_node(Traits::UP).is_null());
  assert(left_top_back.adjacent_node(Traits::BACK).is_null());

  std::cout << octree.root()[Traits::LEFT_BOTTOM_BACK] << std::endl;

  auto right_top_back_of_left_bottom_back = octree.root()[Traits::LEFT_BOTTOM_BACK][Traits::RIGHT_TOP_BACK];
  assert(octree.root()[Traits::LEFT_BOTTOM_BACK][Traits::LEFT_TOP_BACK] == right_top_back_of_left_bottom_back.adjacent_node(Traits::LEFT));
  assert(octree.root()[Traits::RIGHT_BOTTOM_BACK] == right_top_back_of_left_bottom_back.adjacent_node(Traits::RIGHT));
  assert(!right_top_back_of_left_bottom_back.adjacent_node(Traits::RIGHT).is_null());
  assert(!right_top_back_of_left_bottom_back.adjacent_node(Traits::UP).is_null());
  assert(!right_top_back_of_left_bottom_back.adjacent_node(Traits::DOWN).is_null());
  assert(right_top_back_of_left_bottom_back.adjacent_node(Traits::BACK).is_null());
  assert(!right_top_back_of_left_bottom_back.adjacent_node(Traits::FRONT).is_null());

  return 0;
}
