#include <vector>
#include <cassert>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection.h>

using Kernel   = CGAL::Simple_cartesian<double>;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

int main() {

  const std::vector< std::pair<Point_3, Vector_3> > points = {
    std::make_pair(Point_3(0.1, 0.0, 0.0), Vector_3(0.1, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 0.0, 0.0), Vector_3(0.5, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 0.0, 0.0), Vector_3(0.9, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 1.0, 0.0), Vector_3(0.1, 1.0, 1.0)),
    std::make_pair(Point_3(0.5, 1.0, 0.0), Vector_3(0.5, 1.0, 1.0)),
    std::make_pair(Point_3(0.9, 1.0, 0.0), Vector_3(0.9, 1.0, 1.0)),
    std::make_pair(Point_3(0.1, 2.0, 0.0), Vector_3(0.1, 2.0, 1.0)),
    std::make_pair(Point_3(0.5, 2.0, 0.0), Vector_3(0.5, 2.0, 1.0)),
    std::make_pair(Point_3(0.9, 2.0, 0.0), Vector_3(0.9, 2.0, 1.0))
  };

  std::cout << "* number of input points: " << points.size() << std::endl;
  assert(points.size() == 9);

  std::vector< std::vector<std::size_t> > regions;
  CGAL::Shape_detection::region_growing_planes(
    points, std::back_inserter(regions), CGAL::parameters::neighbor_radius(1.5));
  std::cout << "* number of found planar regions: " << regions.size() << std::endl;
  assert(regions.size() == 1);

  regions.clear();
  // CGAL::Shape_detection::region_growing_lines(
  //   points, std::back_inserter(regions), CGAL::parameters::neighbor_radius(0.5));
  std::cout << "* number of found linear regions: " << regions.size() << std::endl;
  assert(regions.size() == 3);

  return EXIT_SUCCESS;
}
