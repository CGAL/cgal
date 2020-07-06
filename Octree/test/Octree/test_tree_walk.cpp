#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Walker_criterion.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

bool test_preorder_1_node() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto tree_walker = CGAL::Octree::Walker::Preorder();
  auto first = tree_walker.first(&octree.root());
  auto nodes = octree.nodes(first, tree_walker);

  // Check each item in the range
  auto iter = nodes.begin();
  if (!(*iter == octree.root()))
    return false;

  return true;
}

bool test_preorder_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto tree_walker = CGAL::Octree::Walker::Preorder();
  auto first = tree_walker.first(&octree.root());
  auto nodes = octree.nodes(first, tree_walker);


  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());
  for (int i = 0; i < 8; ++i) {
    iter++;
    if (!(*iter == octree.root()[i]))
      return false;
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
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto tree_walker = CGAL::Octree::Walker::Preorder();
  auto first = tree_walker.first(&octree.root());
  auto nodes = octree.nodes(first, tree_walker);

  // Check each item in the range
  auto iter = nodes.begin();
  if(!(*iter == octree.root()))
    return false;
  iter++;
  if(!(*iter == octree.root()[0]))
    return false;
  iter++;
  if(!(*iter == octree.root()[1]))
    return false;
  iter++;
  if(!(*iter == octree.root()[2]))
    return false;
  iter++;
  if(!(*iter == octree.root()[3]))
    return false;
  for (int i = 0; i < 8; ++i) {
    iter++;
    if(!(*iter == octree.root()[3][i]))
      return false;
  }
  iter++;
  if(!(*iter == octree.root()[4]))
    return false;
  iter++;
  if(!(*iter == octree.root()[5]))
    return false;
  iter++;
  if(!(*iter == octree.root()[6]))
    return false;
  iter++;
  if(!(*iter == octree.root()[7]))
    return false;
  for (int i = 0; i < 8; ++i) {
    iter++;
    if(!(*iter == octree.root()[7][i]))
      return false;
  }

  return true;
}

int main(void) {

  test_preorder_1_node();
  test_preorder_9_nodes();
  test_preorder_25_nodes();

  return 0;
}
