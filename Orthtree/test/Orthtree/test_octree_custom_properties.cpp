#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>
#include <CGAL/point_generators_3.h>

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

  Octree tree({points, points.point_map()});

  // Testing built in node properties
  typename Octree::Property_map<typename Octree::Node_data> data_prop = tree.get_node_property<typename Octree::Node_data>("contents");

  auto prop2 = tree.get_or_add_node_property("test", int(0));
  assert(prop2.second);

  auto prop3 = tree.get_node_property_if_exists<int>("test");
  assert(prop3.has_value());

  auto prop4 = tree.get_node_property_if_exists<std::string>("test");
  assert(!prop4.has_value());

  auto prop5 = tree.get_or_add_node_property("test", int(0));
  assert(!prop5.second);

  auto prop6 = tree.add_node_property("test2", std::string());

  // Default value should be respected
  auto &node_int_property = tree.add_node_property<int>("int", 5);
  assert(node_int_property[tree.root()] == 5);

  // Changes to individual nodes should be respected
  node_int_property[tree.root()] = 0;
  assert(node_int_property[tree.root()] == 0);

  // Expanding the tree; new nodes should be assigned the default value
  tree.refine(10, 1);
  for (auto n : tree.traverse<CGAL::Orthtrees::Preorder_traversal<Octree>>()) {
    // Everything but the root will have the default value
    if (!tree.is_root(n)) assert(node_int_property[n] == 5);
  }
  // The root should have preserved its custom value
  assert(node_int_property[tree.root()] == 0);




  return EXIT_SUCCESS;
}