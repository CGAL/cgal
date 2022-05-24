#include "include/utils.h"
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <boost/iterator/function_output_iterator.hpp>

// Typedefs.
using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;

using Input_range  = CGAL::Point_set_3<Point_3>;
using Output_range = CGAL::Point_set_3<Point_3>;
using Point_map    = typename Input_range::Point_map;
using Normal_map   = typename Input_range::Vector_map;

using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query, Region_type>;
using Point_inserter = utils::Insert_point_colored_by_region_index<Input_range, Output_range, Point_map>;

int main(int argc, char *argv[]) {

  // Load xyz data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  std::ifstream in(is_default_input ? CGAL::data_file_path("points_3/building.xyz") : argv[1]);

  CGAL::IO::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  const bool with_normal_map = true;
  Input_range input_range(with_normal_map);
  in >> input_range;
  in.close();
  std::cout << "* number of input points: " << input_range.size() << std::endl;
  assert(is_default_input && input_range.size() == 8075);

  // Default parameter values for the data file building.xyz.
  const std::size_t k               = 12;
  const FT          max_distance    = FT(2);
  const FT          max_angle       = FT(20);
  const std::size_t min_region_size = 50;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    input_range, CGAL::parameters::
    k_neighbors(k).
    point_map(input_range.point_map()));

  Region_type region_type(
    input_range,
    CGAL::parameters::
    maximum_distance(max_distance).
    maximum_angle(max_angle).
    minimum_region_size(min_region_size).
    point_map(input_range.point_map()).
    normal_map(input_range.normal_map()));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  // Run the algorithm.
  Output_range output_range;
  std::size_t number_of_regions = 0;
  Point_inserter inserter(
     input_range, input_range.point_map(),
    output_range, number_of_regions);
  region_growing.detect(
    boost::make_function_output_iterator(inserter));
  std::cout << "* number of found planes: " << number_of_regions << std::endl;
  assert(is_default_input && number_of_regions == 7);

  // Save regions to a file.
  const std::string fullpath = (argc > 2 ? argv[2] : "planes_point_set_3.ply");
  std::ofstream out(fullpath);
  out << output_range;
  out.close();

  // Get all unassigned points.
  std::vector<std::size_t> unassigned_items;
  region_growing.unassigned_items(std::back_inserter(unassigned_items));
  std::cout << "* number of unassigned points: " << unassigned_items.size() << std::endl;
  assert(is_default_input && unassigned_items.size() == 535);

  // Store all unassigned points.
  std::vector<Point_3> unassigned_points;
  unassigned_points.reserve(unassigned_items.size());
  for (const std::size_t index : unassigned_items) {
    const auto& key = *(input_range.begin() + index);
    const Point_3& point = get(input_range.point_map(), key);
    unassigned_points.push_back(point);
  }
  return EXIT_SUCCESS;
}
