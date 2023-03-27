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
  assert(octree.adjacent_node(octree.root(), 0) == nullptr);
  assert(octree.adjacent_node(octree.root(), 1) == nullptr);
  assert(octree.adjacent_node(octree.root(), 2) == nullptr);
  assert(octree.adjacent_node(octree.root(), 3) == nullptr);
  assert(octree.adjacent_node(octree.root(), 4) == nullptr);
  assert(octree.adjacent_node(octree.root(), 5) == nullptr);

  // Left Top Front node should have siblings to the Right, Down, and Back
  auto left_top_back = octree.children(octree.root())[Traits::LEFT_TOP_BACK];

  assert(&octree.children(octree.root())[Traits::RIGHT_TOP_BACK] == octree.adjacent_node(left_top_back, Traits::RIGHT));
  assert(
    &octree.children(octree.root())[Traits::LEFT_BOTTOM_BACK] == octree.adjacent_node(left_top_back, Traits::DOWN));
  assert(&octree.children(octree.root())[Traits::LEFT_TOP_FRONT] == octree.adjacent_node(left_top_back, Traits::FRONT));
  assert(octree.adjacent_node(left_top_back, Traits::LEFT) == nullptr);
  assert(octree.adjacent_node(left_top_back, Traits::UP) == nullptr);
  assert(octree.adjacent_node(left_top_back, Traits::BACK) == nullptr);

  std::cout << octree.children(octree.root())[Traits::LEFT_BOTTOM_BACK] << std::endl;

  auto right_top_back_of_left_bottom_back =
    octree.children(octree.children(octree.root())[Traits::LEFT_BOTTOM_BACK])[Traits::RIGHT_TOP_BACK];
  assert(&octree.children(octree.children(octree.root())[Traits::LEFT_BOTTOM_BACK])[Traits::LEFT_TOP_BACK] ==
         octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::LEFT));
  assert(&octree.children(octree.root())[Traits::RIGHT_BOTTOM_BACK] ==
         octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::RIGHT));
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::RIGHT) != nullptr);
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::UP) != nullptr);
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::DOWN) != nullptr);
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::FRONT) != nullptr);

  // A node at the back of the tree should have no neighbor to its back
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::BACK) == nullptr);

  return 0;
}
