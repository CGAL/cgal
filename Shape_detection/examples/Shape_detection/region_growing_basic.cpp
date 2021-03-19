#include <vector>
#include <cassert>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection.h>

using Kernel  = CGAL::Simple_cartesian<float>;
using Point_3 = typename Kernel::Point_3;

int main() {

  const std::vector<Point_3> points = {
    Point_3(0, 0, 0), Point_3(1, 0, 0), Point_3(2, 0, 0),
    Point_3(3, 0, 0), Point_3(4, 0, 0), Point_3(0, 1, 0),
    Point_3(1, 1, 0), Point_3(2, 1, 0), Point_3(3, 1, 0),
    Point_3(4, 1, 0), Point_3(0, 2, 0), Point_3(1, 2, 0),
    Point_3(2, 2, 0), Point_3(3, 2, 0), Point_3(4, 2, 0)
  };
  std::cout << "* number of input points: " << points.size() << std::endl;
  assert(points.size() == 15);

  std::vector< std::vector<std::size_t> > regions;
  // CGAL::Shape_detection::region_growing_planes(points, std::back_inserter(regions));
  std::cout << "* number of found regions: " << regions.size() << std::endl;
  assert(regions.size() == 1);

  return EXIT_SUCCESS;
}
