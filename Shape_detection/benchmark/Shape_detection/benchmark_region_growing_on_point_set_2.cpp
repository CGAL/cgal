// STL includes.
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

namespace SD = CGAL::Shape_detection;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

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

using Timer  = CGAL::Timer;
using Region = std::vector<std::size_t>;

void benchmark_region_growing_on_point_set_2(
  const std::size_t test_count,
  const Input_range& input_range,
  const FT sphere_radius,
  const FT distance_threshold,
  const FT angle_threshold,
  const std::size_t min_region_size) {

  // Create instances of the parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    sphere_radius);

  Region_type region_type(
    input_range,
    distance_threshold, angle_threshold, min_region_size);

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  // Run the algorithm.
  Timer timer;
  std::vector<Region> regions;

  timer.start();
  region_growing.detect(std::back_inserter(regions));
  timer.stop();

  // Compute the number of points assigned to all found regions.
  std::size_t number_of_assigned_points = 0;
  for (const auto& region : regions)
    number_of_assigned_points += region.size();

  std::vector<std::size_t> unassigned_points;
  region_growing.unassigned_items(std::back_inserter(unassigned_points));

  // Print statistics.
  std::cout << "Test #"                          << test_count                << std::endl;
  std::cout << "  sphere_radius = "              << sphere_radius             << std::endl;
  std::cout << "  min_region_size = "            << min_region_size           << std::endl;
  std::cout << "  distance_threshold = "         << distance_threshold        << std::endl;
  std::cout << "  angle_threshold = "            << angle_threshold           << std::endl;
  std::cout << "  -----"                                                      << std::endl;
  std::cout << "  Time elapsed: "                << timer.time()              << std::endl;
  std::cout << "  Number of detected regions: "  << regions.size()            << std::endl;
  std::cout << "  Number of assigned points: "   << number_of_assigned_points << std::endl;
  std::cout << "  Number of unassigned points: " << unassigned_points.size()  << std::endl;
  std::cout << std::endl << std::endl;
}

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_2.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
  in.close();

  // Default parameter values for the data file point_set_2.xyz.
  const FT          distance_threshold = FT(45) / FT(10);
  const FT          angle_threshold    = FT(45);
  const std::size_t min_region_size    = 5;

  // Run benchmarks.
  benchmark_region_growing_on_point_set_2(1, input_range, FT(1),
  distance_threshold, angle_threshold, min_region_size);

  benchmark_region_growing_on_point_set_2(2, input_range, FT(3),
  distance_threshold, angle_threshold, min_region_size);

  benchmark_region_growing_on_point_set_2(3, input_range, FT(6),
  distance_threshold, angle_threshold, min_region_size);

  benchmark_region_growing_on_point_set_2(4, input_range, FT(9),
  distance_threshold, angle_threshold, min_region_size);

  return EXIT_SUCCESS;
}
