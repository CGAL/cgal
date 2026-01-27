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
#include <CGAL/Shape_detection/Region_growing/Point_set.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_point_set_2(int argc, char *argv[]) {

  using FT       = typename Kernel::FT;
  using Point_2  = typename Kernel::Point_2;
  using Vector_2 = typename Kernel::Vector_2;

  using Point_with_normal = std::pair<Point_2, Vector_2>;
  using Input_range       = std::vector<Point_with_normal>;
  using Item              = typename Input_range::const_iterator;
  using Deref_map         = CGAL::Dereference_property_map<const Point_with_normal, Item>;
  using Point_map         = CGAL::Compose_property_map<Deref_map,
                                                      CGAL::First_of_pair_property_map<Point_with_normal>>;
  using Normal_map        = CGAL::Compose_property_map<Deref_map,
                                                      CGAL::Second_of_pair_property_map<Point_with_normal>>;

  using Neighbor_query = SD::Point_set::Sphere_neighbor_query<Kernel, Item, Point_map>;
  using Region_type    = SD::Point_set::Least_squares_line_fit_region<Kernel, Item, Point_map, Normal_map>;
  using Region_growing = SD::Region_growing<Neighbor_query, Region_type>;

  // Default parameter values.
  const FT          sphere_radius      = FT(5);
  const FT          distance_threshold = FT(45) / FT(10);
  const FT          angle_threshold    = FT(45);
  const std::size_t min_region_size    = 5;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : CGAL::data_file_path("points_3/buildings_outline.xyz"));
  CGAL::IO::set_ascii_mode(in);
  assert(in);

  FT a, b, c, d, e, f;
  Input_range input_range;
  while (in >> a >> b >> c >> d >> e >> f) {
    input_range.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
  }
  in.close();
  assert(input_range.size() == 3634);

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    CGAL::parameters::sphere_radius(sphere_radius));

  Region_type region_type(
    CGAL::parameters::
    maximum_distance(distance_threshold).
    maximum_angle(angle_threshold).
    minimum_region_size(min_region_size));

  // Run region growing.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  std::vector<typename Region_growing::Primitive_and_region> regions;
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == 72);
  for (const auto& region : regions)
    assert(region_type.is_valid_region(region.second));

  std::vector<typename Region_growing::Item> unassigned_points;
  region_growing.unassigned_items(input_range, std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 86);
  return true;
}

int main(int argc, char *argv[]) {

  using SC = CGAL::Simple_cartesian<double>;
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

  // ------>

  bool sc_test_success = true;
  if (!test_region_growing_on_point_set_2<SC>(argc, argv))
    sc_test_success = false;
  std::cout << "rg_points2, sc_test_success: " << sc_test_success << std::endl;
  assert(sc_test_success);

  // ------>

  bool epick_test_success = true;
  if (!test_region_growing_on_point_set_2<EPICK>(argc, argv))
    epick_test_success = false;
  std::cout << "rg_points2, epick_test_success: " << epick_test_success << std::endl;
  assert(epick_test_success);

  // ------>

  bool epeck_test_success = true;
  if (!test_region_growing_on_point_set_2<EPECK>(argc, argv))
    epeck_test_success = false;
  std::cout << "rg_points2, epeck_test_success: " << epeck_test_success << std::endl;
  assert(epeck_test_success);

  const bool success = sc_test_success && epick_test_success && epeck_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
