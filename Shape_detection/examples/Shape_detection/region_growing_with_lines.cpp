#include <vector>
#include <cassert>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

// Typedefs.
using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Point_3  = typename Kernel::Point_3;
using Vector_2 = typename Kernel::Vector_2;
using Vector_3 = typename Kernel::Vector_3;

using Point_with_normal_2 = std::pair<Point_2, Vector_2>;
using Point_set_2         = std::vector<Point_with_normal_2>;
using Point_map_2         = CGAL::First_of_pair_property_map<Point_with_normal_2>;
using Normal_map_2        = CGAL::Second_of_pair_property_map<Point_with_normal_2>;

using Point_with_normal_3 = std::pair<Point_3, Vector_3>;
using Point_set_3         = std::vector<Point_with_normal_3>;
using Point_map_3         = CGAL::First_of_pair_property_map<Point_with_normal_3>;
using Normal_map_3        = CGAL::Second_of_pair_property_map<Point_with_normal_3>;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  FT a, b, c, d, e, f;
  Point_set_2 point_set_2;
  Point_set_3 point_set_3;
  while (in >> a >> b >> c >> d >> e >> f) {
    point_set_2.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
    point_set_3.push_back(
      std::make_pair(Point_3(a, b, c), Vector_3(d, e, f)));
  }
  in.close();
  std::cout << "* number of 2D input points: " << point_set_2.size() << std::endl;
  std::cout << "* number of 3D input points: " << point_set_3.size() << std::endl;
  assert(point_set_2.size() == 3634);
  assert(point_set_3.size() == 3634);

  // Run the algorithm both on 2D and 3D points.
  std::vector< std::vector<std::size_t> > lines_2;
  // CGAL::Shape_detection::region_growing_lines(
  //   point_set_2, std::back_inserter(lines_2)
  //   CGAL::parameters::
  //   neighbor_radius(FT(5)).distance_threshold(FT(45) / FT(10)).
  //   angle_threshold(FT(45)).min_region_size(5),
  //   Point_map_2(), Normal_map_2());
  std::cout << "* number of found 2D lines: " << lines_2.size() << std::endl;
  assert(lines_2.size() == 65);

  std::vector< std::vector<std::size_t> > lines_3;
  // CGAL::Shape_detection::region_growing_lines(
  //   point_set_3, std::back_inserter(lines_3)
  //   CGAL::parameters::
  //   neighbor_radius(FT(5)).distance_threshold(FT(45) / FT(10)).
  //   angle_threshold(FT(45)).min_region_size(5),
  //   Point_map_3(), Normal_map_3());
  std::cout << "* number of found 3D lines: " << lines_3.size() << std::endl;
  assert(lines_3.size() == 65);

  return EXIT_SUCCESS;
}
