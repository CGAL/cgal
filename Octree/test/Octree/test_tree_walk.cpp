#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Tree_walker_criterion.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

int test_preorder() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto tree_walker = CGAL::Preorder();
  auto first = tree_walker.first(&octree.root());
  auto nodes = octree.nodes(first, tree_walker);

  for (auto &n : nodes)
    std::cout << n;

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert(*iter == octree.root()[i]);
  }

  return 0;
}

int main(void) {

  test_preorder();
//  test_preorder_print();
//  test_postorder_print();
//  test_leaves_print();

  return 0;
}
