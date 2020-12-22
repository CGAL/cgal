#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

#include <fstream>
#include <iostream>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection_3.h>

// Type declarations.
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// In Shape_detection_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Shape_detection_traits
<Kernel, Pwn_vector, Point_map, Normal_map>                  Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection_3::Region_growing_depr<Traits> Region_growing;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

// This program both works for RANSAC and Region Growing.
// This example is using deprecated code!
// Please update your code to the new version using other examples!
template<typename ShapeDetection>
int run(const char* filename) {

  // Points with normals.
  Pwn_vector points;

  // Load a point set from a file.
  // read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream(filename);

  if (!stream ||
    !CGAL::read_xyz_points(stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map()))) {

    std::cout << "Error: cannot read the file cube.pwn" << std::endl;
    return EXIT_FAILURE;
  }

  // Instantiate a shape detection engine.
  ShapeDetection shape_detection;

  // Provide the input data.
  shape_detection.set_input(points);

  // Register planar shapes via template method.
  shape_detection.template add_shape_factory<Plane>();

  // Detect registered shapes with default parameters.
  shape_detection.detect();

  // Print number of detected shapes.
  std::cout << shape_detection.shapes().end() - shape_detection.shapes().begin()
  << " shapes detected." << std::endl;

  return EXIT_SUCCESS;
}

int main (int argc, char** argv) {

  if (argc > 1 && std::string(argv[1]) == "-r") {
    std::cout << "Efficient RANSAC" << std::endl;
    return run<Efficient_ransac> ((argc > 2) ? argv[2] : "data/cube.pwn");
  }
  std::cout << "Region Growing" << std::endl;
  return run<Region_growing> ((argc > 1) ? argv[1] : "data/cube.pwn");
}
