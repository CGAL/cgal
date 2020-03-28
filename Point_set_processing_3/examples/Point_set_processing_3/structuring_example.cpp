#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/property_map.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>

#include <CGAL/structure_point_set.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// Efficient RANSAC types
typedef CGAL::Shape_detection::Efficient_RANSAC_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>              Traits;
typedef CGAL::Shape_detection::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection::Plane<Traits>               Plane;

int main (int argc, char** argv)
{
  // Points with normals.
  Pwn_vector points;

  // Loading point set from a file.
  std::ifstream stream(argc>1 ? argv[1] : "data/cube.pwn");

  if (!stream ||
    !CGAL::read_xyz_points(stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map())))
  {
      std::cerr << "Error: cannot read file cube.pwn" << std::endl;
      return EXIT_FAILURE;
  }

  std::cerr << points.size() << " point(s) read." << std::endl;

  // Shape detection
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();
  ransac.detect();

  Efficient_ransac::Plane_range planes = ransac.planes();

  Pwn_vector structured_pts;

  CGAL::structure_point_set (points,
                             planes,
                             std::back_inserter (structured_pts),
                             0.015, // epsilon for structuring points
                             CGAL::parameters::point_map (Point_map()).
                             normal_map (Normal_map()).
                             plane_map (CGAL::Shape_detection::Plane_map<Traits>()).
                             plane_index_map (CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)));

  std::cerr << structured_pts.size ()
            << " structured point(s) generated." << std::endl;

  std::ofstream out ("out.pwn");
  CGAL::write_xyz_points (out, structured_pts,
                          CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()));
  out.close();

  return EXIT_SUCCESS;
}
