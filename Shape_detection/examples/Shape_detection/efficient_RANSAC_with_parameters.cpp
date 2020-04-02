#include <fstream>
#include <iostream>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>

// Type declarations.
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::FT                                           FT;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

typedef CGAL::Shape_detection::Efficient_RANSAC_traits
<Kernel, Pwn_vector, Point_map, Normal_map>             Traits;
typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection::Cone<Traits>             Cone;
typedef CGAL::Shape_detection::Cylinder<Traits>         Cylinder;
typedef CGAL::Shape_detection::Plane<Traits>            Plane;
typedef CGAL::Shape_detection::Sphere<Traits>           Sphere;
typedef CGAL::Shape_detection::Torus<Traits>            Torus;

int main(int argc, char** argv) {

  // Points with normals.
  Pwn_vector points;

  // Load point set from a file.
  std::ifstream stream((argc > 1) ? argv[1] : "data/cube.pwn");

  if (!stream ||
    !CGAL::read_xyz_points(
      stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map()))) {

    std::cerr << "Error: cannot read file cube.pwn!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << points.size() << " points" << std::endl;

  // Instantiate shape detection engine.
  Efficient_ransac ransac;

  // Provide input data.
  ransac.set_input(points);

  // Register shapes for detection.
  ransac.add_shape_factory<Plane>();
  ransac.add_shape_factory<Sphere>();
  ransac.add_shape_factory<Cylinder>();
  ransac.add_shape_factory<Cone>();
  ransac.add_shape_factory<Torus>();

  // Set parameters for shape detection.
  Efficient_ransac::Parameters parameters;

  // Set probability to miss the largest primitive at each iteration.
  parameters.probability = 0.05;

  // Detect shapes with at least 200 points.
  parameters.min_points = 200;

  // Set maximum Euclidean distance between a point and a shape.
  parameters.epsilon = 0.002;

  // Set maximum Euclidean distance between points to be clustered.
  parameters.cluster_epsilon = 0.01;

  // Set maximum normal deviation.
  // 0.9 < dot(surface_normal, point_normal);
  parameters.normal_threshold = 0.9;

  // Detect shapes.
  ransac.detect(parameters);

  // Print number of detected shapes and unassigned points.
  std::cout << ransac.shapes().end() - ransac.shapes().begin()
  << " detected shapes, "
  << ransac.number_of_unassigned_points()
  << " unassigned points." << std::endl;

  // Efficient_ransac::shapes() provides
  // an iterator range to the detected shapes.
  Efficient_ransac::Shape_range shapes = ransac.shapes();
  Efficient_ransac::Shape_range::iterator it = shapes.begin();

  while (it != shapes.end()) {

    // Get specific parameters depending on the detected shape.
    if (Plane* plane = dynamic_cast<Plane*>(it->get())) {

      Kernel::Vector_3 normal = plane->plane_normal();
      std::cout << "Plane with normal " << normal << std::endl;

      // Plane shape can also be converted to the Kernel::Plane_3.
      std::cout << "Kernel::Plane_3: " <<
      static_cast<Kernel::Plane_3>(*plane) << std::endl;

    } else if (Cylinder* cyl = dynamic_cast<Cylinder*>(it->get())) {

      Kernel::Line_3 axis = cyl->axis();
      FT radius = cyl->radius();

      std::cout << "Cylinder with axis "
      << axis << " and radius " << radius << std::endl;

    } else {

      // Print the parameters of the detected shape.
      // This function is available for any type of shape.
      std::cout << (*it)->info() << std::endl;
    }

    // Proceed with the next detected shape.
    it++;
  }
  return EXIT_SUCCESS;
}
