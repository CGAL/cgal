#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::vector<Point> Point_vector;
typedef CGAL::Octree<Kernel, Point_vector> Octree;

int main() {

  // Here, our point set is a vector
  Point_vector points;

  // Add a few points to the vector, most of which are in one region
  points.emplace_back(1, 1, 1);
  points.emplace_back(2, 1, -11);
  points.emplace_back(2, 1, 1);
  points.emplace_back(1, -2, 1);
  points.emplace_back(1, 1, 1);
  points.emplace_back(-1, 1, 1);
  points.emplace_back(-1.1, 1, 1);
  points.emplace_back(-1.01, 1, 1);
  points.emplace_back(-1.001, 1, 1);
  points.emplace_back(-1.0001, 1, 1);
  points.emplace_back(-1.0001, 1, 1);

  // Create an octree from the points
  Octree octree(points);

  // Build the octree with a small bucket size, so we get a deep node
  octree.refine(10, 2);

  // Print out the tree
  std::cout << "\nUn-graded tree" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << octree << std::endl;

  // Grade the tree to eliminate large jumps in depth
  octree.grade();

  // Print out the tree again
  std::cout << "\nGraded tree" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << octree << std::endl;

  return EXIT_SUCCESS;
}
