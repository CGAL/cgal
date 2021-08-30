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

int main(int argc, char *argv[]) {
  bool success = true;

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::IO::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file point_set_2.xyz!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    assert(false);
    return EXIT_FAILURE;
  }

  Input_range input_range;
  FT a, b, c, d, e, f;

  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));

  in.close();

  assert(input_range.size() == 3634);
  if (input_range.size() != 3634)
    success = false;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range);
  Region_type region_type(
    input_range);

  // Run region growing.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  // Test data.
  std::vector<std::size_t> unassigned_points;
  region_growing.unassigned_items(std::back_inserter(unassigned_points));

  assert(unassigned_points.size() == 3634);
  if (unassigned_points.size() != 3634)
    success = false;

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));

  assert(regions.size() != 0);
  if (regions.size() == 0)
    success = false;

  unassigned_points.clear();
  region_growing.unassigned_items(std::back_inserter(unassigned_points));

  assert(unassigned_points.size() != 3634);
  if (unassigned_points.size() == 3634)
    success = false;

  const std::size_t num_regions = regions.size();
  regions.clear();
  region_growing.detect(std::back_inserter(regions));

  assert(regions.size() == num_regions);
  if (regions.size() != num_regions)
    success = false;

  std::cout << "basic_test_success: " << success << std::endl;
  assert(success);

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
