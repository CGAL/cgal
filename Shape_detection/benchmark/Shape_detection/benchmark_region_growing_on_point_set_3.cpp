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
using Point_3  = typename Kernel::Point_3;
using Vector_2 = typename Kernel::Vector_2;
using Vector_3 = typename Kernel::Vector_3;

using Point_with_normal = std::pair<Point_3, Vector_3>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Neighbor_query = SD::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = SD::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type>;

using Timer  = CGAL::Timer;
using Region = std::vector<std::size_t>;

void create_input_range(
  const std::size_t num_copies,
  const Input_range& input,
  Input_range& output,
  const bool save = false) {

  const Point_2 a = Point_2(-0.25147, -0.49995);
  const Point_2 b = Point_2( 0.25147, -0.49995);

  output.reserve(num_copies * input.size());
  for (std::size_t i = 0; i < num_copies; ++i) {

    const FT x1 = i * a.x();
    const FT y1 = i * a.y();

    const FT x2 = i * b.x();
    const FT y2 = i * b.y();

    const Point_2 p1 = Point_2(x1, y1);
    const Point_2 p2 = Point_2(x2, y2);

    const Vector_2 tr = Vector_2(p1, p2);
    for (const auto& item : input) {

      const Point_3& p = item.first;
      const Vector_3& n = item.second;

      const Point_3 q = Point_3(
        p.x() + tr.x(),
        p.y() + tr.y(),
        p.z());

      output.push_back(std::make_pair(q, n));
    }
  }

  if (save) {
    const std::string path =
    "path_to_out_folder/";
    const std::string fullpath =
    path + "bench_point_set_3-" + std::to_string(num_copies) + ".xyz";

    std::ofstream out(fullpath);
    for (const auto& item : output)
      out << item.first << " " << item.second << std::endl;
    out.close();
  }
}

void benchmark_region_growing_on_point_set_3(
  const std::size_t num_copies,
  const Input_range& input,
  const std::size_t k,
  const FT distance_threshold,
  const FT angle_threshold,
  const std::size_t min_region_size) {

  Input_range input_range;
  create_input_range(num_copies, input, input_range);

  // Create instances of the parameter classes.
  Neighbor_query neighbor_query(
    input_range,
    k);

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

  std::cout << "Time ( " << input_range.size() <<  " points ): "
  << timer.time() << " seconds" << std::endl;
}

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_3.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_3.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  Input_range input;
  Point_3 p; Vector_3 n;

  while (in >> p >> n)
    input.push_back(std::make_pair(p, n));
  in.close();

  // Default parameter values for the data file point_set_3.xyz.
  const std::size_t k                  = 12;
  const FT          distance_threshold = FT(2);
  const FT          angle_threshold    = FT(20);
  const std::size_t min_region_size    = 25;

  // Run benchmarks.
  std::cout << std::endl;

  benchmark_region_growing_on_point_set_3(1, input,
  k, distance_threshold, angle_threshold, min_region_size);

  benchmark_region_growing_on_point_set_3(2, input,
  k, distance_threshold, angle_threshold, min_region_size);

  benchmark_region_growing_on_point_set_3(3, input,
  k, distance_threshold, angle_threshold, min_region_size);

  benchmark_region_growing_on_point_set_3(4, input,
  k, distance_threshold, angle_threshold, min_region_size);

  std::cout << std::endl << std::endl;
  return EXIT_SUCCESS;
}
