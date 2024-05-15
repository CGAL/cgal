#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Orthtree/Traversals.h>
#include <CGAL/Point_set_3.h>

#include <CGAL/Simple_cartesian.h>
#include <cassert>


using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;
using Preorder_traversal = CGAL::Orthtrees::Preorder_traversal<Octree>;
using Level_traversal = CGAL::Orthtrees::Level_traversal<Octree>;

bool test_preorder_1_node() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Create the range
  auto nodes = octree.traverse<Preorder_traversal>();

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());

  return true;
}

bool test_preorder_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Create the range
  auto nodes = octree.traverse<Preorder_traversal>();

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert(*iter == octree.node(i));
  }

  return true;
}

bool test_level_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);

  // Create the range
  auto nodes = octree.traverse<Level_traversal>(static_cast<std::size_t>(1));

  // Check each item in the range
  auto iter = nodes.begin();
  for (int i = 0; i < 8; ++i) {
    assert(*iter == octree.node(i));
    iter++;
  }

  return true;
}

bool test_preorder_25_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({1, 1, 1});
  points.insert({1, 1, 2});
  points.insert({1, 1, 3});
  points.insert({1, 1, 4});

  // Create the octree
  Octree octree(points, points.point_map());
  octree.refine(10, 1);
  std::cout << octree << std::endl;

  // Create the range
  auto nodes = octree.traverse<Preorder_traversal>();

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());
  iter++;
  assert(*iter == octree.node(0));
  iter++;
  assert(*iter == octree.node(1));
  iter++;
  assert((*iter == octree.node(2)));
  iter++;
  assert(*iter == octree.node(3));
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert(*iter == octree.node(3, i));
  }
  iter++;
  assert((*iter == octree.node(4)));
  iter++;
  assert((*iter == octree.node(5)));
  iter++;
  assert((*iter == octree.node(6)));
  iter++;
  assert((*iter == octree.node(7)));
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert(*iter == octree.node(7, i));
  }

  return true;
}

int main(void) {

  test_preorder_1_node();
  test_preorder_9_nodes();
  test_level_9_nodes();
  test_preorder_25_nodes();

  return 0;
}
