#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Simple_cartesian.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;
using Traits = Octree::Traits;

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
  assert(!octree.adjacent_node(octree.root(), 0));
  assert(!octree.adjacent_node(octree.root(), 1));
  assert(!octree.adjacent_node(octree.root(), 2));
  assert(!octree.adjacent_node(octree.root(), 3));
  assert(!octree.adjacent_node(octree.root(), 4));
  assert(!octree.adjacent_node(octree.root(), 5));

  // Left Top Front node should have siblings to the Right, Down, and Back
  auto left_top_back = octree.node(Traits::LEFT_TOP_BACK);

  assert(octree.node(Traits::RIGHT_TOP_BACK) ==
         *octree.adjacent_node(left_top_back, Traits::RIGHT));
  assert(octree.node(Traits::LEFT_BOTTOM_BACK) ==
         *octree.adjacent_node(left_top_back, Traits::DOWN));
  assert(octree.node(Traits::LEFT_TOP_FRONT) ==
         *octree.adjacent_node(left_top_back, Traits::FRONT));
  assert(!octree.adjacent_node(left_top_back, Traits::LEFT));
  assert(!octree.adjacent_node(left_top_back, Traits::UP));
  assert(!octree.adjacent_node(left_top_back, Traits::BACK));

  auto right_top_back_of_left_bottom_back = octree.node(Traits::LEFT_BOTTOM_BACK, Traits::RIGHT_TOP_BACK);

  assert(
    octree.node(Traits::LEFT_BOTTOM_BACK, Traits::LEFT_TOP_BACK) ==
    octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::LEFT)
  );
  assert(
    octree.node(Traits::RIGHT_BOTTOM_BACK) ==
    octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::RIGHT)
  );
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::RIGHT).has_value());
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::UP).has_value());
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::DOWN).has_value());
  assert(octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::FRONT).has_value());

  // A node at the back of the tree should have no neighbor to its back
  assert(!octree.adjacent_node(right_top_back_of_left_bottom_back, Traits::BACK));

  return 0;
}
