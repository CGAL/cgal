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
using Point_2 = Kernel::Point_2;
using Vector_2 = Kernel::Vector_2;

using Point_set_3 = CGAL::Point_set_3<Point_3, Vector_3>;
using Point_set_2 = CGAL::Point_set_3<Point_2, Vector_2>;
using Point_map = Point_set_2::Point_map;
using Normal_map = Point_set_2::Vector_map;

namespace Shape_detection = CGAL::Shape_detection::Point_set;

using Neighbor_query = Shape_detection::K_neighbor_query
  <Kernel, Point_set_2, Point_map>;
using Region_type = Shape_detection::Least_squares_circle_fit_region
  <Kernel, Point_set_2, Point_map, Normal_map>;
using Sorting = Shape_detection::Least_squares_circle_fit_sorting
  <Kernel, Point_set_2, Neighbor_query, Point_map>;
using Region_growing = CGAL::Shape_detection::Region_growing
  <Point_set_2, Neighbor_query, Region_type, typename Sorting::Seed_map>;

int main (int argc, char** argv)
{
  std::ifstream ifile (argc > 1 ? argv[1] : CGAL::data_file_path("points_3/circles.ply"));
  Point_set_3 points3;
  ifile >> points3;

  std::cerr << points3.size() << " points read" << std::endl;

  // Input should have normals
  assert (points3.has_normal_map());

  Point_set_2 points;
  points.add_normal_map();
  for (Point_set_3::Index idx : points3)
  {
    const Point_3& p = points3.point(idx);
    const Vector_3& n = points3.normal(idx);
    points.insert (Point_2 (p.x(), p.y()), Vector_2 (n.x(), n.y()));
  }

  // Default parameters for data/circles.ply
  const std::size_t k = 12;
  const double tolerance = 0.01;
  const double max_angle = 10.;
  const std::size_t min_region_size = 20;

  // No constraint on radius
  const double min_radius = 0.;
  const double max_radius = std::numeric_limits<double>::infinity();

  Neighbor_query neighbor_query(points, k, points.point_map());
  Region_type region_type(points, tolerance, max_angle, min_region_size,
                          min_radius, max_radius,
                          points.point_map(), points.normal_map());

  // Sort indices
  Sorting sorting(points, neighbor_query, points.point_map());
  sorting.sort();

  Region_growing region_growing(points, neighbor_query, region_type, sorting.seed_map());

  // Add maps to get colored output
  Point_set_3::Property_map<unsigned char>
    red = points3.add_property_map<unsigned char>("red", 0).first,
    green = points3.add_property_map<unsigned char>("green", 0).first,
    blue = points3.add_property_map<unsigned char>("blue", 0).first;

  CGAL::Random random;

  std::size_t nb_circles = 0;
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
        ++ nb_circles;
      }));
  timer.stop();

  std::cerr << nb_circles << " circles detected in "
            << timer.time() << " seconds" << std::endl;

  // Save in colored_circles.ply
  std::ofstream out ("colored_circles.ply");
  CGAL::IO::set_binary_mode (out);
  out << points3;

  return EXIT_SUCCESS;
}
