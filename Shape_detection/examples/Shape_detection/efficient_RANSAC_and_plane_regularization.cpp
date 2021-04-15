#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <fstream>
#include <iostream>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Regularization/regularize_planes.h>

// Type declarations.
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

typedef CGAL::Shape_detection::Efficient_RANSAC_traits
<Kernel, Pwn_vector, Point_map, Normal_map>             Traits;
typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection::Plane<Traits>            Plane;

int main(int argc, char** argv) {

  Pwn_vector points;
  std::ifstream stream(argc > 1 ? argv[1] : "data/cube.pwn");

  if (!stream ||
    !CGAL::read_xyz_points(
      stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map()))) {

    std::cerr << "Error: cannot read file cube.pwn!" << std::endl;
    return EXIT_FAILURE;
  }

  // Call RANSAC shape detection with planes.
  Efficient_ransac efficient_ransac;
  efficient_ransac.set_input(points);
  efficient_ransac.add_shape_factory<Plane>();
  efficient_ransac.detect();

  Efficient_ransac::Plane_range planes = efficient_ransac.planes();

  // Regularize detected planes.
  CGAL::regularize_planes(points,
                          Point_map(),
                          planes,
                          CGAL::Shape_detection::Plane_map<Traits>(),
                          CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                          true,  // regularize parallelism
                          true,  // regularize orthogonality
                          false, // do not regularize coplanarity
                          true,  // regularize Z-symmetry (default)
                          10);   // 10 degrees of tolerance for parallelism / orthogonality

  return EXIT_SUCCESS;
}
