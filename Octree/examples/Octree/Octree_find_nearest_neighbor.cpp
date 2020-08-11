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
  Octree octree(points, points.point_map());

  // Build the octree
  octree.refine(10, 20);

  // Find the nearest points to a few locations
  std::vector<Point> points_to_find = {
          {0, 0, 0},
          {1, 1, 1},
          {-1, -1, -1},
          {-0.46026, -0.25353, 0.32051},
          {-0.460261, -0.253533, 0.320513}
  };
  for (auto p : points_to_find) {

    // The nearest points will be placed in this vector
    std::vector<Point> nearest_points;

    // k=1 to find the single closest point
    octree.nearest_k_neighbors(p, 1, std::back_inserter(nearest_points));

    std::cout << "the nearest point to (" << p << ") is (" << nearest_points[0] << ")" << std::endl;
  }

  return EXIT_SUCCESS;
}
