#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Timer.h> 

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Plane_shape.h>
#include <CGAL/Cylinder_shape.h>
#include <CGAL/Cone_shape.h>
#include <CGAL/Sphere_shape.h>
#include <CGAL/Torus_shape.h>

#include <iostream>
#include <fstream>

// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::vector<Point_with_normal> Point_list;
typedef CGAL::Identity_property_map<Point_with_normal> Point_pmap;
typedef CGAL::Normal_of_point_with_normal_pmap<Kernel> Normal_pmap;

// In Shape_detection_traits_3 the used types, i.e. Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_traits_3<Kernel,
  Point_list::iterator, Point_pmap, Normal_pmap> ShapeDetectionTraits;
typedef CGAL::Shape_detection_3<ShapeDetectionTraits> Shape_detection;


int main(int argc, char **argv) {
  Point_list points;

  // Loading a point set from file. 
  // read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream("cube.xyz");

  if (!stream ||
    !CGAL::read_xyz_points_and_normals(stream,
    std::back_inserter(points),
    Normal_pmap())) {
      std::cerr << "Error: cannot read file cube.xyz" << std::endl;
      return EXIT_FAILURE;
  }

  // Creation of the instance and providing the input data.
  Shape_detection sd(points.begin(),
    points.end(), Point_pmap(), Normal_pmap());

  // Shapes to be searched for are registered
  // by using the template Shape_factory
  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Plane_shape<ShapeDetectionTraits> >);

  // The actual shape detection with default parameters.
  sd.detect();

  // Print results.
  std::cout << sd.number_of_shapes() << " primitives" << std::endl;

  return EXIT_SUCCESS;
}
