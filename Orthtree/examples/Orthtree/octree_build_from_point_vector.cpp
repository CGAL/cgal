#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::list<Point> Point_vector;

typedef CGAL::Octree<Kernel, Point_vector> Octree;

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
  octree.refine(10, 20);

  // Print out the tree
  std::cout << octree;

  return EXIT_SUCCESS;
}
