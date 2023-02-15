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

using Plane_region = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Plane_sorting = SD::Point_set::Least_squares_plane_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Sphere_region = SD::Point_set::Least_squares_sphere_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Sphere_sorting = SD::Point_set::Least_squares_sphere_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Cylinder_region = SD::Point_set::Least_squares_cylinder_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Cylinder_sorting = SD::Point_set::Least_squares_cylinder_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map, Normal_map>;

template <typename Region_type, typename Sorting, typename SortingCode,
          typename RegionCode, typename AssertionCode>
bool test (int argc, char** argv, const SortingCode& sc, const RegionCode& reg, const AssertionCode& assertion)
{
  using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : CGAL::data_file_path("points_3/point_set_3.xyz"));
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

  const std::size_t k = 12;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    k,
    input_range.point_map());

  Region_type region_type = reg(input_range);

  // Sort indices.
  Sorting sorting = sc(input_range, neighbor_query);
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type,
    sorting.seed_map());

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  region_growing.release_memory();

  bool result = assertion(regions);
  assert (result);
  std::cout << "exact_inexact_test_success: " << result << std::endl;
  return result;
}

int main(int argc, char *argv[]) {

  bool success =
    test<Plane_region, Plane_sorting>
    (argc, argv,
     [](const auto& input_range, auto& neighbor_query) -> Plane_sorting
     {
       return Plane_sorting (input_range, neighbor_query,
                             input_range.point_map());
     },
     [](const auto& input_range) -> Plane_region
     {
       // Default parameter values for the data file point_set_3.xyz.
       const FT          distance_threshold = FT(2);
       const FT          angle_threshold    = FT(20);
       const std::size_t min_region_size    = 50;
       return Plane_region
         (input_range,
          distance_threshold, angle_threshold, min_region_size,
          input_range.point_map(), input_range.normal_map());
     },
     [](const auto& r) -> bool {
       std::cout << "- num regions planes: " << r.size() << std::endl;
       return (r.size() >= 6 && r.size() <= 8);
      });
  if (!success)
    return EXIT_FAILURE;

  success =
    test<Sphere_region, Sphere_sorting>
    (argc, argv,
     [](const auto& input_range, auto& neighbor_query) -> Sphere_sorting
     {
       return Sphere_sorting (input_range, neighbor_query,
                              input_range.point_map());
     },
     [](const auto& input_range) -> Sphere_region
     {
       const double tolerance = 0.01;
       const double max_angle = 10.;
       const std::size_t min_region_size = 50;
       // No constraint on radius
       const double min_radius = 0.;
       const double max_radius = std::numeric_limits<double>::infinity();
       return Sphere_region
         (input_range, tolerance, max_angle, min_region_size,
          min_radius, max_radius,
          input_range.point_map(), input_range.normal_map());
     },
     [](const auto& r) -> bool {
        std::cout << "- num regions spheres: " << r.size() << std::endl;
        return (r.size() > 10 && r.size() < 90);
      });
  if (!success)
    return EXIT_FAILURE;

  success =
    test<Cylinder_region, Cylinder_sorting>
    (argc, argv,
     [](const auto& input_range, auto& neighbor_query) -> Cylinder_sorting
     {
       return Cylinder_sorting (input_range, neighbor_query,
                                input_range.point_map(), input_range.normal_map());
     },
     [](const auto& input_range) -> Cylinder_region
     {
       const double tolerance = 0.05;
       const double max_angle = 5.;
       const std::size_t min_region_size = 200;
       // No constraint on radius
       const double min_radius = 0.;
       const double max_radius = std::numeric_limits<double>::infinity();
       return Cylinder_region
         (input_range, tolerance, max_angle, min_region_size,
          min_radius, max_radius,
          input_range.point_map(), input_range.normal_map());
     },
     [](const auto& r) -> bool {
        std::cout << "- num regions cylinders: " << r.size() << std::endl;
        return (r.size() > 2 && r.size() < 30);
      });
  if (!success)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
