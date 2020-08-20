#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Traversal.h>
#include <CGAL/Simple_cartesian.h>

#include <cassert>
#include <CGAL/point_generators_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef std::vector<Point> Point_vector;
typedef CGAL::Octree::Octree<Point_vector> Octree;

int main(void) {


  return EXIT_SUCCESS;
}
