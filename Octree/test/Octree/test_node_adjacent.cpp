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
  std::cout << left_top_back << std::endl;
  left_top_back.adjacent(0);
  left_top_back.adjacent(1);
  left_top_back.adjacent(2);
  left_top_back.adjacent(3);
  left_top_back.adjacent(4);
  left_top_back.adjacent(5);

  assert(nullptr != left_top_back.adjacent(Node::RIGHT));
  assert(nullptr != left_top_back.adjacent(Node::DOWN));
  assert(nullptr != left_top_back.adjacent(Node::FRONT));
  assert(nullptr == left_top_back.adjacent(Node::LEFT));
  assert(nullptr == left_top_back.adjacent(Node::UP));
  assert(nullptr == left_top_back.adjacent(Node::BACK));

  return 0;
}
