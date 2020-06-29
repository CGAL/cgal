#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Tree_walker_criterion.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <CGAL/Octree/Octree_node_range.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

Point_set create_example_point_collection() {

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

  return points;
}

int test_preorder_print() {

  auto points = create_example_point_collection();

  auto point_map = points.point_map();
  Octree octree(points, point_map);
  octree.refine(10, 1);

  std::cout << octree;

  return 0;
}

int test_postorder_print() {

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  auto points = create_example_point_collection();

  auto point_map = points.point_map();
  Octree octree(points, point_map);
  octree.refine(10, 1);

  auto tree_walker = CGAL::Postorder();
  octree.print(std::cout, tree_walker.first(&octree.root()), tree_walker);

  return 0;
}

int test_leaves_print() {

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

  auto points = create_example_point_collection();

  auto point_map = points.point_map();
  Octree octree(points, point_map);
  octree.refine(10, 1);

  auto tree_walker = CGAL::Leaves();
  octree.print(std::cout, tree_walker.first(&octree.root()), tree_walker);

  return 0;
}

int main(void) {

  test_preorder_print();
  test_postorder_print();
  test_leaves_print();

  return 0;
}
