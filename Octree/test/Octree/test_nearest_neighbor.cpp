
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

void naive_vs_accelerated() {

  // Create a dataset

  // Choose a random point

  // Use the naive algorithm to find the nearest point in the dataset

  // Do the same using the octree

  // Check that they produce the same answer

  // Check that the octree was faster

}

int main(void) {

  naive_vs_accelerated();

  return EXIT_SUCCESS;
}
