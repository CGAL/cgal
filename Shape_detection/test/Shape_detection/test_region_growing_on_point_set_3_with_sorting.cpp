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
#include <CGAL/Shape_detection/Region_growing/internal/free_functions.h>

namespace SD = CGAL::Shape_detection;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Input_range = CGAL::Point_set_3<Point_3>;
using Point_map   = typename Input_range::Point_map;
using Normal_map  = typename Input_range::Vector_map;
using Neighbor_query = SD::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;

using Plane_region = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Plane_sorting = SD::Point_set::Least_squares_plane_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Sphere_region = SD::Point_set::Least_squares_sphere_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Sphere_sorting = SD::Point_set::Least_squares_sphere_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Cylinder_region = SD::Point_set::Least_squares_cylinder_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Cylinder_sorting = SD::Point_set::Least_squares_cylinder_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map, Normal_map>;

template<
typename Region_type,
typename Sorting,
typename Lambda_region,
typename Lambda_assertion>
bool test(
  int argc, char** argv, const std::string name,
  const Lambda_region& lambda_region,
  const Lambda_assertion& lambda_assertion) {

  using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  // Default parameter values.
  const std::size_t k = 12;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : CGAL::data_file_path("points_3/building.xyz"));
  CGAL::IO::set_ascii_mode(in);
  assert(in);

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);
  in >> input_range;
  in.close();
  assert(input_range.size() == 8075);

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range, CGAL::parameters::
    k_neighbors(k).
    point_map(input_range.point_map()));

  // Sort indices.
  Sorting sorting(
    input_range, neighbor_query,
    CGAL::parameters::
    point_map(input_range.point_map()).
    normal_map(input_range.normal_map()));
  sorting.sort();

  // Run region growing.
  Region_type region_type = lambda_region(input_range);
  Region_growing region_growing(
    input_range, neighbor_query, region_type, sorting.seed_map());

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  region_growing.clear();
  const bool result = lambda_assertion(regions);
  assert(result);

  // Default parameter values.
  const FT          max_distance    = FT(2);
  const FT          max_angle       = FT(20);
  const std::size_t min_region_size = 50;

  // Test free functions and stability.
  for (std::size_t k = 0; k < 3; ++k) {
    regions.clear();
    SD::internal::region_growing_planes(
      input_range, std::back_inserter(regions),
      CGAL::parameters::
      maximum_distance(max_distance).
      maximum_angle(max_angle).
      minimum_region_size(min_region_size).
      point_map(input_range.point_map()).
      normal_map(input_range.normal_map()));
    assert(regions.size() == 7);
  }

  const bool success = true;
  std::cout << "rg_" + name + "_sortpoints3, epick_test_success: " << success << std::endl;
  return success;
}

int main(int argc, char *argv[]) {
  bool success = true;

  // Test planes.
  success = test<Plane_region, Plane_sorting>(argc, argv, "planes",
    [](const auto& input_range) -> Plane_region {

       const FT          max_distance    = FT(2);
       const FT          max_angle       = FT(20);
       const std::size_t min_region_size = 50;
       return Plane_region(input_range, CGAL::parameters::
          maximum_distance(max_distance).
          maximum_angle(max_angle).
          minimum_region_size(min_region_size).
          point_map(input_range.point_map()).
          normal_map(input_range.normal_map()));
    },
    [](const auto& region) -> bool {
      std::cout << "- num regions planes: " << region.size() << std::endl;
      return (region.size() >= 6 && region.size() <= 8);
    }
  );

  if (!success) {
    return EXIT_FAILURE;
  }

  // Test spheres.
  success = test<Sphere_region, Sphere_sorting>(argc, argv, "spheres",
    [](const auto& input_range) -> Sphere_region {

       const FT          max_distance    = FT(1) / FT(100);
       const FT          max_angle       = FT(10);
       const std::size_t min_region_size = 50;
       return Sphere_region(input_range, CGAL::parameters::
          maximum_distance(max_distance).
          maximum_angle(max_angle).
          minimum_region_size(min_region_size).
          point_map(input_range.point_map()).
          normal_map(input_range.normal_map()));
    },
    [](const auto& region) -> bool {
      std::cout << "- num regions spheres: " << region.size() << std::endl;
      return (region.size() >= 30 && region.size() <= 77);
    }
  );

  if (!success) {
    return EXIT_FAILURE;
  }

  // Test cylinders.
  success = test<Cylinder_region, Cylinder_sorting>(argc, argv, "cylinders",
    [](const auto& input_range) -> Cylinder_region {

       const FT          max_distance    = FT(1) / FT(20);
       const FT          max_angle       = FT(5);
       const std::size_t min_region_size = 200;
       return Cylinder_region(input_range, CGAL::parameters::
          maximum_distance(max_distance).
          maximum_angle(max_angle).
          minimum_region_size(min_region_size).
          point_map(input_range.point_map()).
          normal_map(input_range.normal_map()));
    },
    [](const auto& region) -> bool {
      std::cout << "- num regions cylinders: " << region.size() << std::endl;
      return (region.size() >= 4 && region.size() <= 14);
    }
  );

  if (!success) {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
