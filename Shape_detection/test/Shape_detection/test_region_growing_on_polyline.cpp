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
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>

namespace SD = CGAL::Shape_detection;

template<class Kernel>
bool test_region_growing_on_polyline(int argc, char *argv[]) {

  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;

  using Polyline_2  = std::vector<Point_2>;
  using Polyline_3  = std::vector<Point_3>;
  using Point_map_2 = CGAL::Identity_property_map<Point_2>;
  using Point_map_3 = CGAL::Identity_property_map<Point_3>;

  using Neighbor_query_2 = SD::Polyline::One_ring_neighbor_query<Kernel, Polyline_2>;
  using Region_type_2    = SD::Polyline::Least_squares_line_fit_region<Kernel, Polyline_2, Point_map_2>;
  using Region_growing_2 = SD::Region_growing<Polyline_2, Neighbor_query_2, Region_type_2>;

  using Neighbor_query_3 = SD::Polyline::One_ring_neighbor_query<Kernel, Polyline_3>;
  using Region_type_3    = SD::Polyline::Least_squares_line_fit_region<Kernel, Polyline_3, Point_map_3>;
  using Region_growing_3 = SD::Region_growing<Polyline_3, Neighbor_query_3, Region_type_3>;

  // Default parameter values.
  const FT distance_threshold = FT(45) / FT(10);
  const FT angle_threshold    = FT(45);

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/polyline_3.polylines.txt");
  CGAL::set_ascii_mode(in);
  assert(in);

  // Create 3D polyline.
  Polyline_3 polyline_3;
  std::size_t n = std::size_t(-1);
  in >> n;
  polyline_3.reserve(n);
  while (n--) {
    Point_3 point; in >> point;
    polyline_3.push_back(point);
    assert(in.good());
  }
  in.close();
  assert(polyline_3.size() == 249);

  // Create 3D parameter classes.
  Neighbor_query_3 neighbor_query_3(polyline_3);
  Region_type_3 region_type_3(
    polyline_3,
    CGAL::parameters::
    max_distance(distance_threshold).
    max_angle(angle_threshold));

  // Run 3D region growing.
  Region_growing_3 region_growing_3(
    polyline_3, neighbor_query_3, region_type_3);

  std::vector< std::vector<std::size_t> > regions;
  region_growing_3.detect(std::back_inserter(regions));
  assert(regions.size() == 12);
  for (const auto& region : regions)
    assert(region_type_3.is_valid_region(region));

  std::vector<std::size_t> unassigned_points;
  region_growing_3.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 0);

  // Create 2D polyline.
  std::vector<std::size_t> indices(polyline_3.size());
  std::iota(indices.begin(), indices.end(), 0);
  const auto plane = SD::internal::create_plane(
    polyline_3, Point_map_3(), indices, Kernel()).first;

  Polyline_2 polyline_2;
  polyline_2.reserve(polyline_3.size());
  for (const auto& point : polyline_3) {
    const auto p3 = plane.projection(point);
    const auto p2 = plane.to_2d(p3);
    polyline_2.push_back(p2);
  }
  assert(polyline_2.size() == polyline_3.size());

  // Create 2D parameter classes.
  Neighbor_query_2 neighbor_query_2(polyline_2);
  Region_type_2 region_type_2(
    polyline_2,
    CGAL::parameters::
    max_distance(distance_threshold).
    max_angle(angle_threshold));

  // Run 2D region growing.
  Region_growing_2 region_growing_2(
    polyline_2, neighbor_query_2, region_type_2);

  regions.clear();
  region_growing_2.detect(std::back_inserter(regions));
  assert(regions.size() == 4);
  for (const auto& region : regions)
    assert(region_type_2.is_valid_region(region));

  unassigned_points.clear();
  region_growing_2.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 0);

  return true;
}

int main(int argc, char *argv[]) {

  using SC = CGAL::Simple_cartesian<double>;
  using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
  using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

  // ------>

  bool sc_test_success = true;
  if (!test_region_growing_on_polyline<SC>(argc, argv))
    sc_test_success = false;
  std::cout << "rg_poly, sc_test_success: " << sc_test_success << std::endl;
  assert(sc_test_success);

  // ------>

  bool epick_test_success = true;
  if (!test_region_growing_on_polyline<EPICK>(argc, argv))
    epick_test_success = false;
  std::cout << "rg_poly, epick_test_success: " << epick_test_success << std::endl;
  assert(epick_test_success);

  // ------>

  bool epeck_test_success = true;
  if (!test_region_growing_on_polyline<EPECK>(argc, argv))
    epeck_test_success = false;
  std::cout << "rg_poly, epeck_test_success: " << epeck_test_success << std::endl;
  assert(epeck_test_success);

  const bool success = sc_test_success && epick_test_success && epeck_test_success;
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
