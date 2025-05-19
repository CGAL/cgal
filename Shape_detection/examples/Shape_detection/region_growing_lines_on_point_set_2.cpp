#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>

#include "include/utils.h"

#include <boost/range/irange.hpp>

// Typedefs.
using Kernel   = CGAL::Simple_cartesian<double>;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Point_set_2       = std::vector<Point_with_normal>;

// we use Compose_property_map as the property maps are expected to operate on the item type
// std::size_t passed as parameter to Sphere_neighbor_query and Least_squares_line_fit_region
using Point_map      = CGAL::Compose_property_map<CGAL::Random_access_property_map<Point_set_2>,
                                                  CGAL::First_of_pair_property_map<Point_with_normal> >;
using Normal_map     = CGAL::Compose_property_map<CGAL::Random_access_property_map<Point_set_2>,
                                                  CGAL::Second_of_pair_property_map<Point_with_normal> >;

using Neighbor_query = CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, std::size_t, Point_map>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_line_fit_region<Kernel, std::size_t, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  std::ifstream in(is_default_input ? CGAL::data_file_path("points_3/buildings_outline.xyz") : argv[1]);

  CGAL::IO::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  FT a, b, c, d, e, f;
  Point_set_2 point_set_2;
  while (in >> a >> b >> c >> d >> e >> f) {
    point_set_2.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
  }
  in.close();
  std::cout << "* number of input points: " << point_set_2.size() << std::endl;
  assert(!is_default_input || point_set_2.size() == 3634);

  // Default parameter values for the data file buildings_outline.xyz.
  const FT          sphere_radius   = FT(5);
  const FT          max_distance    = FT(45) / FT(10);
  const FT          max_angle       = FT(45);
  const std::size_t min_region_size = 5;

  Point_map point_map(CGAL::make_random_access_property_map(point_set_2));
  Normal_map normal_map(CGAL::make_random_access_property_map(point_set_2));

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    boost::irange<std::size_t>(0,point_set_2.size()),
    CGAL::parameters::sphere_radius(sphere_radius)
                     .point_map(point_map));

  Region_type region_type(
    CGAL::parameters::
    maximum_distance(max_distance).
    maximum_angle(max_angle).
    minimum_region_size(min_region_size).
    normal_map(normal_map).
    point_map(point_map));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    boost::irange<std::size_t>(0,point_set_2.size()), neighbor_query, region_type);

  // Run the algorithm.
  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));
  std::cout << "* number of found lines: " << regions.size() << std::endl;
  assert(!is_default_input || regions.size() == 72);

  // Save regions to a file.
  const std::string fullpath = (argc > 2 ? argv[2] : "lines_point_set_2.ply");
  utils::save_point_regions_2<Kernel, std::vector<typename Region_growing::Primitive_and_region>, Point_map>(
    regions, fullpath, point_map);
  return EXIT_SUCCESS;
}
