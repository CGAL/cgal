#include "include/utils.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

// Typedefs.
using Kernel   = CGAL::Simple_cartesian<double>;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Neighbor_query = CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_line_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query, Region_type>;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  FT a, b, c, d, e, f;
  Input_range input_range;
  while (in >> a >> b >> c >> d >> e >> f) {
    input_range.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
  }
  in.close();
  std::cout << "* number of input points: " << input_range.size() << std::endl;
  assert(input_range.size() == 3634);

  // Default parameter values for the data file point_set_2.xyz.
  const FT          search_sphere_radius = FT(5);
  const FT          max_distance_to_line = FT(45) / FT(10);
  const FT          max_accepted_angle   = FT(45);
  const std::size_t min_region_size      = 5;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    input_range, CGAL::parameters::neighbor_radius(search_sphere_radius));

  Region_type region_type(
    input_range,
    CGAL::parameters::
    distance_threshold(max_distance_to_line).
    angle_deg_threshold(max_accepted_angle).
    min_region_size(min_region_size));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  std::cout << "* number of found regions: " << regions.size() << std::endl;
  assert(regions.size() == 65);

  // Save regions to a file.
  const std::string fullpath = (argc > 2 ? argv[2] : "regions_point_set_2.ply");
  utils::save_point_regions<Kernel, Input_range, Point_map>(
    input_range, regions, fullpath);
  return EXIT_SUCCESS;
}
