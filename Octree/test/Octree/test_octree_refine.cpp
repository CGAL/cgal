
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;


int test_single_point() {



  return 0;
}

int main(void) {

  int failures = 0;

  failures += test_single_point();

  if (0 == failures)
    return EXIT_SUCCESS;
  return failures;
}
