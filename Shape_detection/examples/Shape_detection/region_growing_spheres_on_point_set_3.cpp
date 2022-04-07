#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <fstream>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Point_set = CGAL::Point_set_3<Point_3>;
using Point_map = typename Point_set::Point_map;
using Normal_map = typename Point_set::Vector_map;

namespace Shape_detection = CGAL::Shape_detection::Point_set;

using Neighbor_query = Shape_detection::K_neighbor_query
  <Kernel, Point_set, Point_map>;
using Region_type = Shape_detection::Least_squares_sphere_fit_region
  <Kernel, Point_set, Point_map, Normal_map>;
using Region_growing = CGAL::Shape_detection::Region_growing
  <Point_set, Neighbor_query, Region_type>;

int main (int argc, char** argv)
{
  std::ifstream ifile (argc > 1 ? argv[1] : CGAL::data_file_path("points_3/spheres.ply"));
  Point_set points;
  ifile >> points;

  std::cerr << points.size() << " points read" << std::endl;

  // Input should have normals
  assert (points.has_normal_map());

  // Default parameters for data/spheres.ply
  const std::size_t k = 12;
  const double tolerance = 0.01;
  const double max_angle = 10.;
  const std::size_t min_region_size = 50;

  // No constraint on radius
  const double min_radius = 0.;
  const double max_radius = std::numeric_limits<double>::infinity();

  Neighbor_query neighbor_query(points, k, points.point_map());
  Region_type region_type(points, tolerance, max_angle, min_region_size,
                          min_radius, max_radius,
                          points.point_map(), points.normal_map());
  Region_growing region_growing(points, neighbor_query, region_type);

  // Add maps to get colored output
  Point_set::Property_map<unsigned char>
    red = points.add_property_map<unsigned char>("red", 0).first,
    green = points.add_property_map<unsigned char>("green", 0).first,
    blue = points.add_property_map<unsigned char>("blue", 0).first;

  CGAL::Random random;

  std::size_t nb_spheres = 0;
  CGAL::Real_timer timer;
  timer.start();
  region_growing.detect
    (boost::make_function_output_iterator
     ([&](const std::vector<std::size_t>& region)
      {
        // Assign a random color to each region
        unsigned char r = static_cast<unsigned char>(random.get_int(64, 192));
        unsigned char g = static_cast<unsigned char>(random.get_int(64, 192));
        unsigned char b = static_cast<unsigned char>(random.get_int(64, 192));
        for (const std::size_t& idx : region)
        {
          red[idx] = r;
          green[idx] = g;
          blue[idx] = b;
        }
        ++ nb_spheres;
      }));
  timer.stop();

  std::cerr << nb_spheres << " spheres detected in "
            << timer.time() << " seconds" << std::endl;

  // Save in colored_spheres.ply
  std::ofstream out ("colored_spheres.ply");
  CGAL::IO::set_binary_mode (out);
  out << points;

  return EXIT_SUCCESS;
}
