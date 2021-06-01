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

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::IO::set_ascii_mode(in);
  assert(in);

  FT a, b, c, d, e, f;
  Input_range input_range;
  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(std::make_pair(Point_2(a, b), Vector_2(d, e)));
  in.close();
  assert(input_range.size() == 3634);

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range, CGAL::parameters::all_default());
  Region_type region_type(
    input_range, CGAL::parameters::all_default());

  // Run region growing.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  std::vector<std::size_t> unassigned_points;
  region_growing.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 3634);

  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  const std::size_t num_regions = regions.size();
  assert(num_regions != 0);

  unassigned_points.clear();
  region_growing.unassigned_items(std::back_inserter(unassigned_points));
  const std::size_t num_unassigned_points = unassigned_points.size();
  assert(num_unassigned_points != 3634);

  regions.clear();
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == num_regions);

  unassigned_points.clear();
  region_growing.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == num_unassigned_points);

  region_growing.clear();

  unassigned_points.clear();
  region_growing.unassigned_items(std::back_inserter(unassigned_points));
  assert(unassigned_points.size() == 3634);

  std::cout << "rg_basic, epick_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}
