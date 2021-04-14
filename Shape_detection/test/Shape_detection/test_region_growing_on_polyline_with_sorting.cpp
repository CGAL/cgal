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
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>
#include <CGAL/Shape_detection/Region_growing/internal/free_functions.h>

namespace SD = CGAL::Shape_detection;
using Kernel = CGAL::Simple_cartesian<double>;

using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using Polyline_2  = std::vector<Point_2>;
using Polyline_3  = std::vector<Point_3>;
using Point_map_2 = CGAL::Identity_property_map<Point_2>;
using Point_map_3 = CGAL::Identity_property_map<Point_3>;

using Neighbor_query_2 = SD::Polyline::One_ring_neighbor_query<Kernel, Polyline_2>;
using Region_type_2    = SD::Polyline::Least_squares_line_fit_region<Kernel, Polyline_2, Point_map_2>;
using Sorting_2        = SD::Polyline::Least_squares_line_fit_sorting<Kernel, Polyline_2, Neighbor_query_2, Point_map_2>;
using Region_growing_2 = SD::Region_growing<Polyline_2, Neighbor_query_2, Region_type_2, typename Sorting_2::Seed_map>;

using Neighbor_query_3 = SD::Polyline::One_ring_neighbor_query<Kernel, Polyline_3>;
using Region_type_3    = SD::Polyline::Least_squares_line_fit_region<Kernel, Polyline_3, Point_map_3>;
using Sorting_3        = SD::Polyline::Least_squares_line_fit_sorting<Kernel, Polyline_3, Neighbor_query_3, Point_map_3>;
using Region_growing_3 = SD::Region_growing<Polyline_3, Neighbor_query_3, Region_type_3, typename Sorting_3::Seed_map>;

int main(int argc, char *argv[]) {

  // Default parameter values.
  const FT max_distance_to_line = FT(45) / FT(10);
  const FT max_accepted_angle   = FT(45);

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
    distance_threshold(max_distance_to_line).
    angle_threshold(max_accepted_angle));

  // Sort indices.
  Sorting_3 sorting_3(
    polyline_3, neighbor_query_3, CGAL::parameters::all_default());
  sorting_3.sort();

  // Run 3D region growing.
  Region_growing_3 region_growing_3(
    polyline_3, neighbor_query_3, region_type_3, sorting_3.seed_map());

  std::vector< std::vector<std::size_t> > regions;
  region_growing_3.detect(std::back_inserter(regions));
  assert(regions.size() == 15);
  for (const auto& region : regions)
    assert(region_type_3.is_valid_region(region));

  std::vector<std::size_t> unassigned_points;
  region_growing_3.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 0);

  // Test free functions and stability.
  for (std::size_t k = 0; k < 3; ++k) {
    regions.clear();
    SD::internal::region_growing_polylines(
      polyline_3, std::back_inserter(regions),
      CGAL::parameters::
      distance_threshold(max_distance_to_line).
      angle_threshold(max_accepted_angle));
    assert(regions.size() == 15);
  }

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
    distance_threshold(max_distance_to_line).
    angle_threshold(max_accepted_angle));

  // Sort indices.
  Sorting_2 sorting_2(
    polyline_2, neighbor_query_2, CGAL::parameters::all_default());
  sorting_2.sort();

  // Run 2D region growing.
  Region_growing_2 region_growing_2(
    polyline_2, neighbor_query_2, region_type_2, sorting_2.seed_map());

  regions.clear();
  region_growing_2.detect(std::back_inserter(regions));
  assert(regions.size() == 5);
  for (const auto& region : regions)
    assert(region_type_2.is_valid_region(region));

  unassigned_points.clear();
  region_growing_2.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 0);

  // Test free functions and stability.
  for (std::size_t k = 0; k < 3; ++k) {
    regions.clear();
    SD::internal::region_growing_polylines(
      polyline_2, std::back_inserter(regions),
      CGAL::parameters::
      distance_threshold(max_distance_to_line).
      angle_threshold(max_accepted_angle));
    assert(regions.size() == 5);
  }

  std::cout << "rg_sortpoly, sc_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}
