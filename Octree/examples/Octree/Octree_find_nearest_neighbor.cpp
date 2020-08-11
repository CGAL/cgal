#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;
typedef CGAL::Octree::Octree<Point_set, Point_map> Octree;

int main(int argc, char **argv) {

  // Point set will be used to hold our points
  Point_set points;

  // Load points from a file.
  std::ifstream stream((argc > 1) ? argv[1] : "data/cube.pwn");
  stream >> points;
  if (0 == points.number_of_points()) {

    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "loaded " << points.number_of_points() << " points" << std::endl;

  // Create an octree from the points
  Point_map point_map = points.point_map();
  Octree octree(points, point_map);

  // Build the octree
  octree.refine(10, 20);

  // Find the nearest points to a few locations
  std::vector<Point> points_to_find = {
          {0, 0, 0},
          {1, 1, 1}
  };
  for (auto p : points_to_find) {

    std::cout << "the nearest point to (" << p << ") is (" << std::endl;
  }

  return EXIT_SUCCESS;
}
