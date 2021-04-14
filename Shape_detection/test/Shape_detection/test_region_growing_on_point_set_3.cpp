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
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_point_set_3(int argc, char *argv[]) {

  using FT      = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  using Input_range = CGAL::Point_set_3<Point_3>;
  using Point_map   = typename Input_range::Point_map;
  using Normal_map  = typename Input_range::Vector_map;

  using Neighbor_query = SD::Point_set::Sphere_neighbor_query<Kernel, Input_range, Point_map>;
  using Region_type    = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
  using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type>;

  // Default parameter values.
  const FT          sphere_radius      = FT(5) / FT(100);
  const FT          distance_threshold = FT(2);
  const FT          angle_threshold    = FT(20);
  const std::size_t min_region_size    = 50;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_3.xyz");
  CGAL::set_ascii_mode(in);
  assert(in);

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);
  in >> input_range;
  in.close();
  assert(input_range.size() == 8075);

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range, CGAL::parameters::
    neighbor_radius(sphere_radius).
    point_map(input_range.point_map()));

  Region_type region_type(
    input_range,
    CGAL::parameters::
    distance_threshold(distance_threshold).
    angle_threshold(angle_threshold).
    min_region_size(min_region_size).
    point_map(input_range.point_map()).
    normal_map(input_range.normal_map()));

  // Run region growing.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == 9);
  for (const auto& region : regions)
    assert(region_type.is_valid_region(region));

  std::vector<std::size_t> unassigned_points;
  region_growing.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 196);
  return true;
}

int main(int argc, char *argv[]) {

  using SC = CGAL::Simple_cartesian<double>;
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;

  // ------>

  bool sc_test_success = true;
  if (!test_region_growing_on_point_set_3<SC>(argc, argv))
    sc_test_success = false;
  std::cout << "rg_points3, sc_test_success: " << sc_test_success << std::endl;
  assert(sc_test_success);

  // ------>

  bool epick_test_success = true;
  if (!test_region_growing_on_point_set_3<EPICK>(argc, argv))
    epick_test_success = false;
  std::cout << "rg_points3, epick_test_success: " << epick_test_success << std::endl;
  assert(epick_test_success);

  // ------>

  bool epeck_test_success = true; // turn it off in case it is slow
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
  if (!test_region_growing_on_point_set_3<EPECK>(argc, argv))
    epeck_test_success = false;
  std::cout << "rg_points3, epeck_test_success: " << epeck_test_success << std::endl;
  assert(epeck_test_success);

  const bool success = sc_test_success && epick_test_success && epeck_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
