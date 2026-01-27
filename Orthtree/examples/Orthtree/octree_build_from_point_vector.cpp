#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>

// Type Declarations
using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Point_vector = std::vector<Point>;
using Octree = CGAL::Octree<Kernel, Point_vector>;

int main() {

  // Here, our point set is a vector
  Point_vector points;

  // Add a few points to the vector
  points.emplace_back(1, 1, 1);
  points.emplace_back(2, 1, -11);
  points.emplace_back(2, 1, 1);
  points.emplace_back(1, -2, 1);
  points.emplace_back(1, 1, 1);
  points.emplace_back(-1, 1, 1);

  // Create an octree from the points
  Octree octree(points);

  // Build the octree
  octree.refine(10, 1);

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
