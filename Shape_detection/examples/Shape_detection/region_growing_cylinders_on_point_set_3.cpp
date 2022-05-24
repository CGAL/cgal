#include "include/utils.h"
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
#include <boost/iterator/function_output_iterator.hpp>

// Typedefs.
using Kernel   = CGAL::Simple_cartesian<double>;
using FT       = Kernel::FT;
using Point_3  = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Point_set  = CGAL::Point_set_3<Point_3>;
using Point_map  = typename Point_set::Point_map;
using Normal_map = typename Point_set::Vector_map;

using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query<Kernel, Point_set, Point_map>;
using Region_type    = CGAL::Shape_detection::Point_set::Least_squares_cylinder_fit_region<Kernel, Point_set, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing<Point_set, Neighbor_query, Region_type>;

int main(int argc, char** argv) {

  // Load ply data either from a local folder or a user-provided file.
  const bool is_default_input = argc > 1 ? false : true;
  std::ifstream in(is_default_input ? CGAL::data_file_path("points_3/cube.pwn") : argv[1]);

  CGAL::IO::set_ascii_mode(in);
  if (!in) {
    std::cerr << "ERROR: cannot read the input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Point_set point_set;
  in >> point_set;
  in.close();
  std::cout << "* number of input points: " << point_set.size() << std::endl;
  assert(is_default_input && point_set.size() == 50000);
  assert(point_set.has_normal_map()); // input should have normals

  // Default parameter values for the data file cuble.pwn.
  const std::size_t k               = 24;
  const FT          max_distance    = FT(1) / FT(20);
  const FT          max_angle       = FT(5);
  const std::size_t min_region_size = 200;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    point_set, CGAL::parameters::k_neighbors(k).point_map(point_set.point_map()));

  Region_type region_type(
    point_set,
    CGAL::parameters::
    maximum_distance(max_distance).
    maximum_angle(max_angle).
    minimum_region_size(min_region_size).
    point_map(point_set.point_map()).
    normal_map(point_set.normal_map()));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    point_set, neighbor_query, region_type);

  // Add maps to get a colored output.
  Point_set::Property_map<unsigned char>
    red   = point_set.add_property_map<unsigned char>("red"  , 0).first,
    green = point_set.add_property_map<unsigned char>("green", 0).first,
    blue  = point_set.add_property_map<unsigned char>("blue" , 0).first;

  // Run the algorithm.
  CGAL::Random random;
  std::size_t num_cylinders = 0;
  region_growing.detect(
    boost::make_function_output_iterator(
      [&](const std::vector<std::size_t>& region) {

        // Assign a random color to each region.
        const unsigned char r = static_cast<unsigned char>(random.get_int(64, 192));
        const unsigned char g = static_cast<unsigned char>(random.get_int(64, 192));
        const unsigned char b = static_cast<unsigned char>(random.get_int(64, 192));
        for (const std::size_t idx : region) {
          red[idx]   = r;
          green[idx] = g;
          blue[idx]  = b;
        }
        ++num_cylinders;
      }
    )
  );
  std::cout << "* number of found cylinders: " << num_cylinders << std::endl;
  assert(is_default_input && num_cylinders == 62);

  // Save regions to a file.
  std::ofstream out("cylinders_point_set_3.ply");
  CGAL::IO::set_ascii_mode(out);
  out << point_set;
  return EXIT_SUCCESS;
}
