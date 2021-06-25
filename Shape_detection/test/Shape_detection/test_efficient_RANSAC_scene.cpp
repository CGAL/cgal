#include "include/test_efficient_RANSAC_generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Regularization/regularize_planes.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>

template <class K>
bool test_scene(int argc, char** argv) {

  typedef typename K::FT                                      FT;
  typedef CGAL::Point_with_normal_3<K>                        Pwn;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef CGAL::Identity_property_map<Pwn>                    Point_map;
  typedef CGAL::Normal_of_point_with_normal_map<K>            Normal_map;

  typedef CGAL::Shape_detection::Efficient_RANSAC_traits<K, Pwn_vector, Point_map, Normal_map> Traits;
  typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;

  typedef typename Efficient_ransac::Point_index_range Point_index_range;

  typedef CGAL::Shape_detection::Plane<Traits>              Plane;
  typedef CGAL::Shape_detection::Cone<Traits>               Cone;
  typedef CGAL::Shape_detection::Cylinder<Traits>           Cylinder;
  typedef CGAL::Shape_detection::Sphere<Traits>             Sphere;
  typedef CGAL::Shape_detection::Torus<Traits>              Torus;

  Pwn_vector points;

  // Load point set from a file.
  // read_points takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  const char* filename = (argc > 1) ? argv[1] : "data/cube.pwn";

  if (!CGAL::IO::read_points(filename, std::back_inserter(points),
                         CGAL::parameters::point_map(Point_map())
                                          .normal_map(Normal_map())))
  {
    std::cerr << "Error: cannot read file cube.pwn" << std::endl;
    return EXIT_FAILURE;
  }

  Efficient_ransac ransac;

  ransac.template add_shape_factory<Cone>();

  ransac.clear_shape_factories();

  ransac.template add_shape_factory<Cone>();
  ransac.template add_shape_factory<Cylinder>();
  ransac.template add_shape_factory<Plane>();
  ransac.template add_shape_factory<Sphere>();
  ransac.template add_shape_factory<Torus>();

  ransac.set_input(points);

  ransac.preprocess();

  ransac.clear_octrees();

  if (!ransac.detect()) {
    std::cout << " aborted" << std::endl;
    return false;
  }

  typename Efficient_ransac::Shape_range shapes = ransac.shapes();

  typename Efficient_ransac::Shape_range::iterator it = shapes.begin();

  FT average_distance = 0;

  // Iterate through all shapes and access each point.
  while (it != shapes.end()) {
    boost::shared_ptr<typename Efficient_ransac::Shape> shape = *it;

    // Sum distances of points to detected shapes.
    FT sum_distances = 0;

    // Iterate through point indices assigned to each detected shape.
    std::vector<std::size_t>::const_iterator
      index_it = (*it)->indices_of_assigned_points().begin();

    while (index_it != (*it)->indices_of_assigned_points().end()) {

      // Retrieve point.
      const Pwn &p = *(points.begin() + (*index_it));

      // Add Euclidean distance between point and shape.
      sum_distances += CGAL::sqrt((*it)->squared_distance(p));

      // Proceed with next point.
      index_it++;
    }

    // Compute average distance.
    average_distance += sum_distances / shape->indices_of_assigned_points().size();

    // Proceed with next detected shape.
    it++;
  }

  // Check coverage. For this scene it should not fall below 75%.
  double coverage = double(points.size() - ransac.number_of_unassigned_points()) / double(points.size());
  if (coverage < 0.75) {
    std::cout << " failed (coverage = " << coverage << " < 0.75)" << std::endl;

    return false;
  }

  // Check average distance. It should not lie above 0.02.
  average_distance = average_distance / shapes.size();
  std::cout << average_distance << " " << std::endl;
  if (average_distance > 0.02) {
    std::cout << " failed" << std::endl;

    return false;
  }

  // Test regularization
  typename Efficient_ransac::Plane_range planes = ransac.planes();
  CGAL::regularize_planes (points,
                           Point_map(),
                           planes,
                           CGAL::Shape_detection::Plane_map<Traits>(),
                           CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes),
                           true, true, true, true,
                           (FT)50., (FT)0.01);

  Point_index_range pts = ransac.indices_of_unassigned_points();

  std::cout << " succeeded" << std::endl;

  return true;
}

int main(int argc, char** argv) {
  bool success = true;

  std::cout << "test_scene<CGAL::Simple_cartesian<float>> ";
  if (!test_scene<CGAL::Simple_cartesian<float> >(argc, argv))
    success = false;

//  std::cout << "test_scene<CGAL::Simple_cartesian<double>> ";
//  if (!test_scene<CGAL::Simple_cartesian<double> >(argc, argv))
//    success = false;
//
//  std::cout << "test_scene<CGAL::Exact_predicates_inexact_constructions_kernel> ";
//  if (!test_scene<CGAL::Exact_predicates_inexact_constructions_kernel>(argc, argv))
//    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
