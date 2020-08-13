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
  auto &left_top_back = octree.root()[Node::LEFT_TOP_BACK];

  assert(octree.root()[Node::RIGHT_TOP_BACK] == *left_top_back.adjacent(Node::RIGHT));
  assert(octree.root()[Node::LEFT_BOTTOM_BACK] == *left_top_back.adjacent(Node::DOWN));
  assert(octree.root()[Node::LEFT_TOP_FRONT] == *left_top_back.adjacent(Node::FRONT));
  assert(nullptr == left_top_back.adjacent(Node::LEFT));
  assert(nullptr == left_top_back.adjacent(Node::UP));
  assert(nullptr == left_top_back.adjacent(Node::BACK));

  std::cout << octree.root()[Node::LEFT_BOTTOM_BACK] << std::endl;

  auto &right_top_back_of_left_bottom_back = octree.root()[Node::LEFT_BOTTOM_BACK][Node::RIGHT_TOP_BACK];
  assert(octree.root()[Node::LEFT_BOTTOM_BACK][Node::LEFT_TOP_BACK] == *right_top_back_of_left_bottom_back.adjacent(Node::LEFT));
  assert(octree.root()[Node::RIGHT_BOTTOM_BACK][Node::LEFT_TOP_BACK] == *right_top_back_of_left_bottom_back.adjacent(Node::RIGHT));
  assert(nullptr != right_top_back_of_left_bottom_back.adjacent(Node::RIGHT));
  assert(nullptr != right_top_back_of_left_bottom_back.adjacent(Node::UP));
  assert(nullptr != right_top_back_of_left_bottom_back.adjacent(Node::DOWN));
  assert(nullptr == right_top_back_of_left_bottom_back.adjacent(Node::BACK));
  assert(nullptr != right_top_back_of_left_bottom_back.adjacent(Node::FRONT));

  return 0;
}
