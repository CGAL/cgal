#include "include/utils.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polyline.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;
using Plane_3 = typename Kernel::Plane_3;

using Polyline_2  = std::vector<Point_2>;
using Polyline_3  = std::vector<Point_3>;
using Point_map_2 = CGAL::Identity_property_map<Point_2>;
using Point_map_3 = CGAL::Identity_property_map<Point_3>;

using Neighbor_query = CGAL::Shape_detection::Polyline::One_ring_neighbor_query<Kernel, Polyline_3>;
using Region_type    = CGAL::Shape_detection::Polyline::Least_squares_line_fit_region<Kernel, Polyline_3, Point_map_3>;
using Region_growing = CGAL::Shape_detection::Region_growing<Polyline_3, Neighbor_query, Region_type>;

int main(int argc, char *argv[]) {

  // Load polyline data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  std::ifstream in(is_default_input ? "data/polyline_3.polylines.txt" : argv[1]);
  CGAL::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Polyline_3 polyline_3;
  std::size_t n = std::size_t(-1);
  in >> n;
  polyline_3.reserve(n);
  while (n--) {
    Point_3 point; in >> point;
    polyline_3.push_back(point);
    if (!in.good()) {
      std::cout << "ERROR: cannot load a polyline!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  in.close();
  std::cout << "* number of input vertices: " << polyline_3.size() << std::endl;
  assert(is_default_input && polyline_3.size() == 249);

  // Default parameter values for the data file polyline_3.polylines.txt.
  const FT max_distance_to_line = FT(45) / FT(10);
  const FT max_accepted_angle   = FT(45);

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(polyline_3);

  Region_type region_type(
    polyline_3,
    CGAL::parameters::
    distance_threshold(max_distance_to_line).
    angle_threshold(max_accepted_angle));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    polyline_3, neighbor_query, region_type);

  // Run the algorithm on a 3D polyline.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  std::cout << "* number of found 3D regions: " << regions.size() << std::endl;
  assert(is_default_input && regions.size() == 12);

  // Save 3D regions to a file.
  // std::string fullpath = (argc > 2 ? argv[2] : "regions_polyline_3.ply");
  // utils::save_point_regions_3<Kernel, Polyline_3, Point_map_3>(
  //   polyline_3, regions, fullpath);

  // Create the 2D polyline.
  Plane_3 plane; Point_3 centroid;
  CGAL::linear_least_squares_fitting_3(
    polyline_3.begin(), polyline_3.end(), plane, centroid,
    CGAL::Dimension_tag<0>(), Kernel(),
    CGAL::Eigen_diagonalize_traits<FT, 3>());

  Polyline_2 polyline_2;
  polyline_2.reserve(polyline_3.size());
  for (const auto& point : polyline_3) {
    const auto p3 = plane.projection(point);
    const auto p2 = plane.to_2d(p3);
    polyline_2.push_back(p2);
  }
  assert(is_default_input && polyline_2.size() == polyline_3.size());

  // Use a free function to get the 2D regions.
  regions.clear();
  CGAL::Shape_detection::internal::region_growing_polylines(
    polyline_2, std::back_inserter(regions),
    CGAL::parameters::
    distance_threshold(max_distance_to_line).
    angle_threshold(max_accepted_angle));
  std::cout << "* number of found 2D regions: " << regions.size() << std::endl;
  assert(is_default_input && regions.size() == 5);

  // Save 2D regions to a file.
  // fullpath = (argc > 2 ? argv[2] : "regions_polyline_2.ply");
  // utils::save_point_regions_2<Kernel, Polyline_2, Point_map_2>(
  //   polyline_2, regions, fullpath);

  return EXIT_SUCCESS;
}
