#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::vector<Point> Point_set;
typedef CGAL::Identity_property_map<Point> Point_map;
typedef CGAL::Octree::Octree<Point_set, Point_map> Octree;

int main(int argc, char **argv) {

  // Here, our point set is a vector
  Point_set points;

  // Add a few points to the vector
  points.emplace_back(1, 1, 1);
  points.emplace_back(2, 1, -11);
  points.emplace_back(2, 1, 1);
  points.emplace_back(1, -2, 1);
  points.emplace_back(1, 1, 1);
  points.emplace_back(-1, 1, 1);

  // Create an octree from the points
  Point_map point_map = Point_map ();
  Octree octree(points, point_map);

  // Build the octree with a small bucket size
  octree.refine(10, 2);

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
