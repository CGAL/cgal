#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

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
typedef Kernel::Point_3 Point;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::vector<Point_with_normal> Point_list;
typedef CGAL::Identity_property_map<Point_with_normal> Point_pmap;
typedef CGAL::Normal_of_point_with_normal_pmap<Kernel> Normal_pmap;


int main(int argc, char **argv) {
  Point_list points;

  // Loading a point set from file. read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream("cube.xyz");
  if (!stream || !CGAL::read_xyz_points_and_normals(stream, std::back_inserter(points), Normal_pmap())) {
    std::cerr << "Error: cannot read file cube.xyz" << std::endl;
    return EXIT_FAILURE;
  }
    
  // For time measurement we take the time before setting up the shape detection.
  clock_t start, end;
  start = clock();
    
  // In Shape_detection_traits_3 the used types, i.e. Point and Vector types as well as iterator type and property maps, are defined.
  typedef CGAL::Shape_detection_traits_3<Kernel, Point_list::iterator, Point_pmap, Normal_pmap> ShapeDetectionTraits;
  typedef CGAL::Shape_detection_3<ShapeDetectionTraits> Shape_detection;

  // Creation of the instance and providing the input data.
  Shape_detection sd(points.begin(), points.end(), Point_pmap(), Normal_pmap());
    
  // Shapes to be searched for are registered by using the template Shape_factory
  sd.add_shape_factory(new CGAL::Shape_factory<CGAL::Plane_shape<ShapeDetectionTraits> >);
  sd.add_shape_factory(new CGAL::Shape_factory<CGAL::Cylinder_shape<ShapeDetectionTraits> >);
    
  // Parameterization of the shape detection using the Parameters structure.
  Shape_detection::Parameters op;
  op.probability = 0.01f;       // 1% probability to miss the largest primitive on each iteration.
  op.min_points = 500;          // Only extract shapes with at least 500 points.
  op.epsilon = 0.002f;          // 0.002 maximum euclidean distance between point and shape.
  op.cluster_epsilon = 0.01f;    // 0.01 maximum euclidean distance between points to be clustered.
  op.normal_threshold = 0.9f;   // 0.9 < dot(surface_normal, point_normal); maximum normal deviation.
  
  // The actual shape detection.
  sd.detect(op);
    
  // Take the time afterwards.
  end = clock();
  std::cout << "time: " << (double(end - start) / CLOCKS_PER_SEC) * 1000 << "ms" << std::endl;
  std::cout << sd.number_of_shapes() << " primitives, " << sd.number_of_unassigned_points() << " unassigned points" << std::endl;
  
  // Shape_detection_3::shapes_begin() provides an iterator to the detected shapes.
  auto it = sd.shapes_begin();
  while (it != sd.shapes_end()) {
    // Using Shape_base::info() for printing the parameters of the detected shape.
    std::cout << (*it)->info() << std::endl;
    it++;
  }

  return EXIT_SUCCESS;
}
