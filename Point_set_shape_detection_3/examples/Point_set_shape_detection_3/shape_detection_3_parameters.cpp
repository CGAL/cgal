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
typedef Kernel::FT                                          FT;
typedef CGAL::Point_with_normal_3<Kernel>                   Point_with_normal;
typedef std::vector<Point_with_normal>                      Pwn_vector;
typedef CGAL::Identity_property_map<Point_with_normal>      Point_pmap;
typedef CGAL::Normal_of_point_with_normal_pmap<Kernel>      Normal_pmap;

// In Shape_detection_traits_3 the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_traits_3<Kernel,
  Pwn_vector::iterator, Point_pmap, Normal_pmap>            Traits;
typedef CGAL::Shape_detection_3<Traits>                     Shape_detection;


int main() 
{
  // Points with normals.
  Pwn_vector points;

  // Loads point set from a file. 
  // read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream("cube.pwn");

  if (!stream ||
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Normal_pmap())) 
  {
    std::cerr << "Error: cannot read file cube.pwn" << std::endl;
    return EXIT_FAILURE;
  }

  // Instantiates shape detection engine and provides input data.
  Shape_detection sd(points.begin(),
    points.end(), Point_pmap(), Normal_pmap());
    
  // Register shapes for detection
	
  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Plane_shape<Traits> >);
		
  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Cylinder_shape<Traits> >);

  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Sphere_shape<Traits> >);
 
  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Cone_shape<Traits> >);
 
  sd.add_shape_factory(new 
    CGAL::Shape_factory<CGAL::Torus_shape<Traits> >);
		

  // Sets parameters for shape detection.
  Shape_detection::Parameters parameters;

  // Sets probability to miss the largest primitive at each iteration.
  parameters.probability = 0.05;
 
  // Detect shapes with at least 500 points.
  parameters.min_points = 500;

  // Sets maximum Euclidean distance between a point and a shape.
  parameters.epsilon = 0.002;
 
  // Sets maximum Euclidean distance between points to be clustered.
  parameters.cluster_epsilon = 0.01;
 
  // Sets maximum normal deviation.
  // 0.9 < dot(surface_normal, point_normal); 
  parameters.normal_threshold = 0.9;   
  
  // Detects shapes
  sd.detect(parameters);

  // Prints number of detected shapes and unassigned points.
   std::cout << sd.number_of_shapes() << " detected shapes, "
     << sd.number_of_unassigned_points()
     << " unassigned points." << std::endl;
  
  // Shape_detection_3::shapes_begin() provides
  // an iterator to the detected shapes.
  Shape_detection::Shape_iterator it = sd.shapes_begin();
  while (it != sd.shapes_end()) {
    const Shape_detection::Shape *shape = *it;
    // Prints the parameters of the detected shape.
    std::cout << (*it)->info() << std::endl;

    // Proceeds with next detected shape.
    it++;
  }

  return EXIT_SUCCESS;
}
