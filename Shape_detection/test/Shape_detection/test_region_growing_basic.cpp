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
#include <CGAL/Shape_detection/Region_growing/internal/free_functions.h>

namespace SD = CGAL::Shape_detection;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Point_3  = typename Kernel::Point_3;
using Vector_2 = typename Kernel::Vector_2;
using Vector_3 = typename Kernel::Vector_3;

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
  CGAL::set_ascii_mode(in);
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
  assert(unassigned_points.size() != 3634);

  regions.clear();
  region_growing.detect(std::back_inserter(regions));
  assert(regions.size() == num_regions);

  const std::vector< std::pair<Point_3, Vector_3> > points_with_normals = {
    std::make_pair(Point_3(0.1, 0.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 0.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 0.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 1.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 1.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 1.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.1, 2.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.5, 2.0, 0.0), Vector_3(0.0, 0.0, 1.0)),
    std::make_pair(Point_3(0.9, 2.0, 0.0), Vector_3(0.0, 0.0, 1.0))
  };
  assert(points_with_normals.size() == 9);

  regions.clear();
  CGAL::Shape_detection::internal::region_growing_planes(
    points_with_normals, std::back_inserter(regions));
  assert(regions.size() == 1);
  assert(regions[0].size() == 9);

  std::cout << "rg_basic, epick_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}
