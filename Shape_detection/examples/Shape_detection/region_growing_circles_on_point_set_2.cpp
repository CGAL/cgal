#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>
#include <boost/iterator/function_output_iterator.hpp>

#include "include/utils.h"

// Typedefs.
using Kernel   = CGAL::Simple_cartesian<double>;
using FT       = Kernel::FT;
using Point_2  = Kernel::Point_2;
using Point_3  = Kernel::Point_3;
using Vector_2 = Kernel::Vector_2;
using Vector_3 = Kernel::Vector_3;

using Point_set_2 = CGAL::Point_set_3<Point_2, Vector_2>;
using Point_set_3 = CGAL::Point_set_3<Point_3, Vector_3>;

using Point_map  = Point_set_2::Point_map;
using Normal_map = Point_set_2::Vector_map;

using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_set_2>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_circle_fit_region_for_point_set<Point_set_2>;
using Sorting        = CGAL::Shape_detection::Point_set::Least_squares_circle_fit_sorting_for_point_set<Point_set_2, Neighbor_query>;
using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

int main(int argc, char** argv) {

  // Load ply data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  std::ifstream in(is_default_input ? CGAL::data_file_path("points_3/circles.ply") : argv[1]);

  CGAL::IO::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Point_set_3 point_set_3;
  in >> point_set_3;
  in.close();
  std::cout << "* number of input points: " << point_set_3.size() << std::endl;
  assert(!is_default_input || point_set_3.size() == 1101);
  assert(point_set_3.has_normal_map()); // input should have normals

  // Create a 2D point set.
  Point_set_2 point_set_2;
  point_set_2.add_normal_map();
  for (const auto& idx : point_set_3) {
    const Point_3& point = point_set_3.point(idx);
    const Vector_3& normal = point_set_3.normal(idx);
    point_set_2.insert(
      Point_2(point.x(), point.y()),
      Vector_2(normal.x(), normal.y()));
  }

  // Default parameter values for the data file circles.ply.
  const std::size_t k               = 16;
  const FT          max_distance    = FT(1) / FT(100);
  const FT          max_angle       = FT(10);
  const std::size_t min_region_size = 20;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(
    point_set_2, CGAL::parameters::k_neighbors(k));

  Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_circle_fit_region(
    point_set_2,
    CGAL::parameters::
    maximum_distance(max_distance).
    maximum_angle(max_angle).
    minimum_region_size(min_region_size));

  // Sort indices.
  Sorting sorting = CGAL::Shape_detection::Point_set::make_least_squares_circle_fit_sorting(
    point_set_2, neighbor_query, CGAL::parameters::k_neighbors(k));
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    point_set_2, sorting.ordered(), neighbor_query, region_type);

  // Add maps to get a colored output.
  Point_set_2::Property_map<unsigned char>
    red   = point_set_2.add_property_map<unsigned char>("red"  , 0).first,
    green = point_set_2.add_property_map<unsigned char>("green", 0).first,
    blue  = point_set_2.add_property_map<unsigned char>("blue" , 0).first;

  // Run the algorithm.
  CGAL::Random random;
  std::size_t num_circles = 0;
  region_growing.detect(
    boost::make_function_output_iterator(
      [&](const std::pair< Region_type::Primitive, typename Region_growing::Region >& region) {

        // Assign a random color to each region.
        const unsigned char r = static_cast<unsigned char>(random.get_int(64, 192));
        const unsigned char g = static_cast<unsigned char>(random.get_int(64, 192));
        const unsigned char b = static_cast<unsigned char>(random.get_int(64, 192));
        for (auto item : region.second) {
          red[item] = r;
          green[item] = g;
          blue[item] = b;
        }
        ++num_circles;
      }
    )
  );
  std::cout << "* number of found circles: " << num_circles << std::endl;
  assert(!is_default_input || num_circles == 10);

  // Save regions to a file.
  std::ofstream out("circles_point_set_2.ply");
  CGAL::IO::set_ascii_mode(out);
  out << point_set_2;
  return EXIT_SUCCESS;
}
