// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cassert>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_point_set_2(int argc, char *argv[]) {

  using FT       = typename Kernel::FT;
  using Point_2  = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;

  using Point_with_normal = std::pair<Point_2, Vector_2>;
  using Input_range       = std::vector<Point_with_normal>;
  using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
  using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

  using Neighbor_query = SD::Point_set::Sphere_neighbor_query<Kernel, Input_range, Point_map>;
  using Region_type    = SD::Point_set::Least_squares_line_fit_region<Kernel, Input_range, Point_map, Normal_map>;
  using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type>;

  // Default parameter values for the data file point_set_2.xyz.
  const FT          sphere_radius      = FT(5);
  const FT          distance_threshold = FT(45) / FT(10);
  const FT          angle_threshold    = FT(45);
  const std::size_t min_region_size    = 5;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_2.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    return false;
  }

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));

  in.close();

  assert(input_range.size() == 3634);
  if (input_range.size() != 3634)
    return false;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    sphere_radius);

  Region_type region_type(
    input_range,
    distance_threshold, angle_threshold, min_region_size);

  // Run region growing.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  // Test data.
  assert(regions.size() >= 63 && regions.size() <= 67);
  if (regions.size() < 63 || regions.size() > 67)
    return false;

  for (const auto& region : regions)
    if (!region_type.is_valid_region(region))
      return false;

  std::vector<std::size_t> unassigned_points;
  region_growing.unassigned_items(std::back_inserter(unassigned_points));

  assert(unassigned_points.size() >= 77 && unassigned_points.size() <= 97);
  if (unassigned_points.size() < 77 || unassigned_points.size() > 97)
    return false;

  return true;
}

int main(int argc, char *argv[]) {

  // ------>

  bool cartesian_double_test_success = true;
  if (!test_region_growing_on_point_set_2< CGAL::Simple_cartesian<double> >(argc, argv))
    cartesian_double_test_success = false;

  std::cout << "cartesian_double_test_success: " << cartesian_double_test_success << std::endl;
  assert(cartesian_double_test_success);

  // ------>

  bool exact_inexact_test_success = true;
  if (!test_region_growing_on_point_set_2<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv))
    exact_inexact_test_success = false;

  std::cout << "exact_inexact_test_success: " << exact_inexact_test_success << std::endl;
  assert(exact_inexact_test_success);

  // ------>

  bool exact_exact_test_success = true;
  if (!test_region_growing_on_point_set_2<CGAL::Exact_predicates_exact_constructions_kernel>(argc, argv))
    exact_exact_test_success = false;

  std::cout << "exact_exact_test_success: " << exact_exact_test_success << std::endl;
  assert(exact_exact_test_success);

  const bool success = cartesian_double_test_success && exact_inexact_test_success && exact_exact_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
