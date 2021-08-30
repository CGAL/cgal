// STL includes.
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_point_set_3(int argc, char *argv[]) {

  using FT       = typename Kernel::FT;
  using Point_3  = typename Kernel::Point_3;

  using Input_range = CGAL::Point_set_3<Point_3>;
  using Point_map   = typename Input_range::Point_map;
  using Normal_map  = typename Input_range::Vector_map;

  using Neighbor_query = SD::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
  using Region_type    = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
  using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type>;

  // Default parameter values for the data file point_set_3.xyz.
  const std::size_t k                  = 12;
  const FT          distance_threshold = FT(2);
  const FT          angle_threshold    = FT(20);
  const std::size_t min_region_size    = 50;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_3.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_3.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    return false;
  }

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);

  in >> input_range;
  in.close();

  assert(input_range.size() == 8075);
  if (input_range.size() != 8075)
    return false;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    k,
    input_range.point_map());

  Region_type region_type(
    input_range,
    distance_threshold, angle_threshold, min_region_size,
    input_range.point_map(), input_range.normal_map());

  // Run region growing.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  // Test data.
  assert(regions.size() >= 6 && regions.size() <= 8);
  if (regions.size() < 6 || regions.size() > 8)
    return false;

  for (const auto& region : regions)
    if (!region_type.is_valid_region(region))
      return false;

  std::vector<std::size_t> unassigned_points;
  region_growing.unassigned_items(std::back_inserter(unassigned_points));

  assert(unassigned_points.size() >= 528 && unassigned_points.size() <= 548);
  if (unassigned_points.size() < 528 || unassigned_points.size() > 548)
    return false;

  return true;
}

int main(int argc, char *argv[]) {

  // ------>

  bool cartesian_double_test_success = true;
  if (!test_region_growing_on_point_set_3< CGAL::Simple_cartesian<double> >(argc, argv))
    cartesian_double_test_success = false;

  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
  assert(cartesian_double_test_success);

  // ------>

  bool exact_inexact_test_success = true;
  if (!test_region_growing_on_point_set_3<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv))
    exact_inexact_test_success = false;

  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
  assert(exact_inexact_test_success);

  const bool success = cartesian_double_test_success && exact_inexact_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
