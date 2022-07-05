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
#include <CGAL/Shape_detection/Region_growing/Polyline.h>
#include <CGAL/Shape_detection/Region_growing/free_functions.h>

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
using Region_growing_2 = SD::Region_growing<Polyline_2, Neighbor_query_2, Region_type_2>;

using Neighbor_query_3 = SD::Polyline::One_ring_neighbor_query<Kernel, Polyline_3>;
using Region_type_3    = SD::Polyline::Least_squares_line_fit_region<Kernel, Polyline_3, Point_map_3>;
using Sorting_3        = SD::Polyline::Least_squares_line_fit_sorting<Kernel, Polyline_3, Neighbor_query_3, Point_map_3>;
using Region_growing_3 = SD::Region_growing<Polyline_3, Neighbor_query_3, Region_type_3>;

int main(int argc, char *argv[]) {

  // Default parameter values.
  const FT distance_threshold = FT(45) / FT(10);
  const FT angle_threshold    = FT(45);

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : CGAL::data_file_path("polylines_3/wavy_circle.polylines.txt"));
  CGAL::IO::set_ascii_mode(in);
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
    maximum_distance(distance_threshold).
    maximum_angle(angle_threshold));

  // Sort indices.
  Sorting_3 sorting_3(polyline_3, neighbor_query_3);
  sorting_3.sort();

  // Run 3D region growing.
  Region_growing_3 region_growing_3(
    polyline_3, neighbor_query_3, region_type_3, sorting_3.ordered());

  Region_growing_3::Result_type regions3;
  region_growing_3.detect(std::back_inserter(regions3));
  assert(regions3.size() == 16);

  Region_growing_3::Unassigned_type unassigned_points;
  region_growing_3.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 1);

  // Test free functions and stability.
  for (std::size_t k = 0; k < 3; ++k) {
    regions3.clear();
    SD::region_growing_polylines(
      polyline_3, std::back_inserter(regions3),
      CGAL::parameters::
      maximum_distance(distance_threshold).
      maximum_angle(angle_threshold));
    assert(regions3.size() == 16);
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
    maximum_distance(distance_threshold).
    maximum_angle(angle_threshold));

  // Sort indices.
  Sorting_2 sorting_2(polyline_2, neighbor_query_2);
  sorting_2.sort();

  // Run 2D region growing.
  Region_growing_2 region_growing_2(
    polyline_2, neighbor_query_2, region_type_2, sorting_2.ordered());

  typename Region_growing_2::Result_type regions2;
  region_growing_2.detect(std::back_inserter(regions2));
  assert(regions2.size() == 5);
  for (const auto& region : regions2)
    assert(region_type_2.is_valid_region(region.second));

  Region_growing_2::Region_map map = region_growing_2.region_map();

  for (std::size_t i = 0; i < regions2.size(); i++)
    for (auto& item : regions2[i].second) {
      if (i != get(map, CGAL::Shape_detection::internal::conditional_deref<Region_growing_2::Item, typename Region_growing_2::Region_map::key_type>()(item))) {
        std::cout << "Region map incorrect" << std::endl;
        assert(false);
      }
    }

  Region_growing_2::Unassigned_type unassigned;
  region_growing_2.unassigned_items(std::back_inserter(unassigned));

  for (auto& item : unassigned) {
    if (std::size_t(-1) != get(map, CGAL::Shape_detection::internal::conditional_deref<Region_growing_2::Item, typename Region_growing_2::Region_map::key_type>()(item))) {
      std::cout << "Region map for unassigned incorrect" << std::endl;
      assert(false);
    }
  }

  Region_growing_2::Unassigned_type unassigned2;
  region_growing_2.unassigned_items(std::back_inserter(unassigned2));
  assert(unassigned2.size() == 0);

  // Test free functions and stability.
  for (std::size_t k = 0; k < 3; ++k) {
    regions2.clear();
    SD::region_growing_polylines(
      polyline_2, std::back_inserter(regions2),
      CGAL::parameters::
      maximum_distance(distance_threshold).
      maximum_angle(angle_threshold));
    assert(regions2.size() == 5);
  }

  std::cout << "rg_sortpoly, sc_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}