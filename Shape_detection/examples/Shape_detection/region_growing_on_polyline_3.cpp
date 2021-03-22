#include "include/utils.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Input_range = std::vector<Point_3>;
using Point_map   = CGAL::Identity_property_map<Point_3>;

using Neighbor_query = CGAL::Shape_detection::Polyline::One_ring_neighbor_query<Kernel, Input_range, Point_map>;
using Region_type    = CGAL::Shape_detection::Polyline::Least_squares_line_fit_region<Kernel, Input_range, Point_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Input_range, Neighbor_query, Region_type>;

int main(int argc, char *argv[]) {

  // Load polyline data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/polyline_3.polylines.txt");
  CGAL::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  std::size_t n = std::size_t(-1);
  in >> n;
  Input_range input_range;
  input_range.reserve(n);
  while (n--) {
    Point_3 point; in >> point;
    input_range.push_back(point);
    if (!in.good()) {
      std::cout << "ERROR: cannot load a polyline!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  in.close();
  std::cout << "* number of input vertices: " << input_range.size() << std::endl;
  assert(input_range.size() == 249);

  // Default parameter values for the data file polyline_3.polylines.txt.
  const FT          max_distance_to_line = FT(45) / FT(10);
  const FT          max_accepted_angle   = FT(45);
  const std::size_t min_region_size      = 5;
  const bool        sort_regions         = true;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(input_range);

  Region_type region_type(
    input_range,
    CGAL::parameters::
    distance_threshold(max_distance_to_line).
    angle_deg_threshold(max_accepted_angle).
    min_region_size(min_region_size).
    sort_regions(sort_regions));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    input_range, neighbor_query, region_type);

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  std::cout << "* number of found regions: " << regions.size() << std::endl;
  assert(regions.size() == 9);

  // Sort the output.
  if (sort_regions) {
    std::sort(regions.begin(), regions.end(),
    [&](const std::vector<std::size_t>& r1, const std::vector<std::size_t>& r2) {
      // Since all indices within each region are sorted, see the region_type,
      // and our input is a polyline, it is enough to sort them by the first index only.
      assert(r1.size() >= min_region_size);
      assert(r2.size() >= min_region_size);
      return r1.front() < r2.front();
    });
    for (std::size_t i = 0; i < regions.size(); ++i) {
      std::cout << "region " << i << ": ";
      for (const std::size_t index : regions[i]) {
        std::cout << index << " ";
      }
      std::cout <<std::endl;
    }
  }

  // Save regions to a file.
  const std::string fullpath = (argc > 2 ? argv[2] : "regions_polyline_3.ply");
  utils::save_point_regions<Kernel, Input_range, Point_map>(
    input_range, regions, fullpath);
  return EXIT_SUCCESS;
}
