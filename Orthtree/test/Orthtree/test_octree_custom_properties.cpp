#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Simple_cartesian.h>

#include <iostream>
#include <cassert>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using FT = Kernel::FT;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;

int main(void) {

  std::size_t nb_pts = 100;
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(nb_pts);
  for (std::size_t i = 0; i < nb_pts; ++i)
    points.insert(*(generator++));

  Octree tree(points, points.point_map());

  // Testing built in node properties
  typename Octree::Property_map<typename Octree::Node_data> data_prop = *tree.property<typename Octree::Node_data>("contents");
  CGAL_USE(data_prop);

  // list of properties
  assert(tree.properties().size() == 5);
  auto prop1 = tree.add_property("test", int(5));
  assert(prop1.second);
  // One property added
  assert(tree.properties().size() == 6);
  // Default value should be respected
  assert(prop1.first[tree.root()] == 5);
  // Changes to individual nodes should be respected
  prop1.first[tree.root()] = 0;
  assert(prop1.first[tree.root()] == 0);

  auto prop2 = tree.add_property("test", int(0));
  assert(!prop2.second);
  assert(tree.properties().size() == 6);

  auto prop3 = tree.add_property("test2", std::string());
  assert(prop3.second);
  assert(tree.properties().size() == 7);

  auto prop4 = tree.property<int>("test");
  assert(prop4.has_value());

  auto prop5 = tree.property<std::string>("test");
  assert(!prop5.has_value());

  // Expanding the tree; new nodes should be assigned the default value
  tree.refine(10, 1);
  for (auto n : tree.traverse<CGAL::Orthtrees::Preorder_traversal<Octree>>()) {
    // Everything but the root will have the default value
    if (!tree.is_root(n)) assert(prop1.first[n] == 5);
  }
  // The root should have preserved its custom value
  assert(prop1.first[tree.root()] == 0);

  tree.remove_property(prop1.first);
  assert(tree.properties().size() == 6);

  return EXIT_SUCCESS;
}
