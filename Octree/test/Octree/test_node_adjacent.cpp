#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Traversal.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;
typedef Octree::Node Node;

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
  assert(nullptr == octree.root().adjacent(0));
  assert(nullptr == octree.root().adjacent(1));
  assert(nullptr == octree.root().adjacent(2));
  assert(nullptr == octree.root().adjacent(3));
  assert(nullptr == octree.root().adjacent(4));
  assert(nullptr == octree.root().adjacent(5));

  // Left Top Front node should have siblings to the Right, Down, and Back
  auto &left_top_front = octree.root()[Node::LEFT_TOP_FRONT];
  assert(nullptr != left_top_front.adjacent(Node::RIGHT));
  assert(nullptr != left_top_front.adjacent(Node::DOWN));
  std::cout << left_top_front << std::endl;
  std::cout << *left_top_front.adjacent(Node::RIGHT) << std::endl;
  assert(nullptr != left_top_front.adjacent(Node::BACK));
  assert(nullptr == left_top_front.adjacent(Node::LEFT));
  assert(nullptr == left_top_front.adjacent(Node::UP));
  assert(nullptr == left_top_front.adjacent(Node::FRONT));

  return 0;
}
