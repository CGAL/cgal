// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Neighbor_query = SD::Point_set::Sphere_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = SD::Point_set::Least_squares_line_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = SD::Point_set::Least_squares_line_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_2.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    assert(false);
    return EXIT_FAILURE;
  }

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
  in.close();

  // Default parameter values for the data file point_set_2.xyz.
  const FT          sphere_radius      = FT(5);
  const FT          distance_threshold = FT(45) / FT(10);
  const FT          angle_threshold    = FT(45);
  const std::size_t min_region_size    = 5;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    sphere_radius);

  Region_type region_type(
    input_range,
    distance_threshold, angle_threshold, min_region_size);

  // Sort indices.
  Sorting sorting(
    input_range, neighbor_query);
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type,
    sorting.seed_map());

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  region_growing.release_memory();
  assert(regions.size() >= 62 && regions.size() <= 66);

  const bool cartesian_double_test_success = (regions.size() >= 62 && regions.size() <= 66);
  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;

  const bool success = cartesian_double_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
