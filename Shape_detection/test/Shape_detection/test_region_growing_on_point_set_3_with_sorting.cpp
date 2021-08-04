// STL includes.
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

// Type declarations.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;

using Input_range = CGAL::Point_set_3<Point_3>;
using Point_map   = typename Input_range::Point_map;
using Normal_map  = typename Input_range::Vector_map;

using Neighbor_query = SD::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = SD::Point_set::Least_squares_plane_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_3.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_3.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    assert(false);
    return EXIT_FAILURE;
  }

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);

  in >> input_range;
  in.close();

  // Default parameter values for the data file point_set_3.xyz.
  const std::size_t k                  = 12;
  const FT          distance_threshold = FT(2);
  const FT          angle_threshold    = FT(20);
  const std::size_t min_region_size    = 50;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    k,
    input_range.point_map());

  Region_type region_type(
    input_range,
    distance_threshold, angle_threshold, min_region_size,
    input_range.point_map(), input_range.normal_map());

  // Sort indices.
  Sorting sorting(
    input_range, neighbor_query,
    input_range.point_map());
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type,
    sorting.seed_map());

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  region_growing.release_memory();
  assert(regions.size() >= 6 && regions.size() <= 8);

  const bool exact_inexact_test_success = (regions.size() >= 6 && regions.size() <= 8);
  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;

  const bool success = exact_inexact_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
