#include "include/utils.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Shape_regularization/regularize_planes.h>

// Typedefs.
using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using Point_with_normal = std::pair<Point_3, Vector_3>;
using Point_map         = CGAL::First_of_pair_property_map<Point_with_normal>;
using Normal_map        = CGAL::Second_of_pair_property_map<Point_with_normal>;
using Pwn_vector        = std::vector<Point_with_normal>;

using Traits    = CGAL::Shape_detection::Efficient_RANSAC_traits<Kernel, Pwn_vector, Point_map, Normal_map>;
using RANSAC    = CGAL::Shape_detection::Efficient_RANSAC<Traits>;
using Plane     = CGAL::Shape_detection::Plane<Traits>;
using Plane_map = CGAL::Shape_detection::Plane_map<Traits>;

int main(int argc, char** argv) {

  // If we want to load a different file, we load it from a path.
  std::string path = CGAL::data_file_path("points_3/cube.pwn");
  if (argc > 1) path = argv[1];

  Pwn_vector points;
  std::ifstream file(path.c_str(), std::ios_base::in);
  CGAL::IO::set_ascii_mode(file);
  file.precision(20);

  if (!file ||
    !CGAL::IO::read_XYZ(
      file,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map()))) {

    std::cerr << "Error: cannot read the file cube.pwn!" << std::endl;
    return EXIT_FAILURE;
  }
  file.close();

  // Call RANSAC shape detection with planes.
  RANSAC efficient_ransac;
  efficient_ransac.set_input(points);
  efficient_ransac.add_shape_factory<Plane>();
  efficient_ransac.detect();

  auto planes = efficient_ransac.planes();

  // Regularize detected planes.
  CGAL::Shape_regularization::Planes::regularize_planes(
    planes,
    points,
    CGAL::parameters::
    plane_map(Plane_map()).
    point_map(Point_map()).
    plane_index_map(
      CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)).
    regularize_coplanarity(false). // do not regularize coplanarity
    maximum_angle(FT(10))); // 10 degrees of tolerance for parallelism / orthogonality

  std::cout << "* all detected planes are regularized" << std::endl;
  return EXIT_SUCCESS;
}
