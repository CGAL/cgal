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
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <CGAL/Shape_detection/Region_growing/internal/free_functions.h>

namespace SD = CGAL::Shape_detection;
using Kernel = CGAL::Simple_cartesian<double>;

using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Vector_2 = typename Kernel::Vector_2;

using Point_with_normal = std::pair<Point_2, Vector_2>;
using Input_range       = std::vector<Point_with_normal>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;

using Neighbor_query = SD::Point_set::K_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = SD::Point_set::Least_squares_line_fit_region<Kernel, Input_range, Point_map, Normal_map>;
using Sorting        = SD::Point_set::Least_squares_line_fit_sorting<Kernel, Input_range, Neighbor_query, Point_map>;
using Region_growing = SD::Region_growing<Input_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  // Load data.
  std::ifstream in(argc > 1 ? argv[1] : "data/point_set_2.xyz");
  CGAL::set_ascii_mode(in);
  assert(in);

  FT a, b, c, d, e, f;
  Input_range input_range;
  while (in >> a >> b >> c >> d >> e >> f)
    input_range.push_back(
      std::make_pair(Point_2(a, b), Vector_2(d, e)));
  in.close();
  assert(input_range.size() == 3634);

  // Default parameter values for the data file point_set_2.xyz.
  const std::size_t k                  = 12;
  const FT          distance_threshold = FT(45) / FT(10);
  const FT          angle_threshold    = FT(45);
  const std::size_t min_region_size    = 5;

  // Create parameter classes.
  Neighbor_query neighbor_query(
    input_range, CGAL::parameters::neighbor_radius(k));

  Region_type region_type(
    input_range,
    CGAL::parameters::
    distance_threshold(distance_threshold).
    angle_threshold(angle_threshold).
    min_region_size(min_region_size));

  // Sort indices.
  Sorting sorting(
    input_range, neighbor_query, CGAL::parameters::all_default());
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type, sorting.seed_map());

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  region_growing.release_memory();
  // std::cout << regions.size() << std::endl;
  assert(regions.size() == 62);

  // Test determenistic behavior and free functions.
  for (std::size_t k = 0; k < 3; ++k) {
    regions.clear();
    SD::internal::region_growing_lines(
      input_range, std::back_inserter(regions),
      CGAL::parameters::
      distance_threshold(distance_threshold).
      angle_threshold(angle_threshold).
      min_region_size(min_region_size));
    assert(regions.size() == 62);
  }

  std::cout << "rg_sortpoints2, sc_test_success: " << true << std::endl;
  return EXIT_SUCCESS;
}
