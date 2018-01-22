#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/Point_set_processing_3.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Shape_detection_3.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Random_points_on_sphere_3<Point> Point_generator;

typedef CGAL::Point_set_3<Point> Point_set;

typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits
<Kernel, Point_set, Point_set::Point_map, Point_set::Vector_map> Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>        Efficient_ransac;
typedef CGAL::Shape_detection_3::Sphere<Traits>                  Sphere;


int main (int, char**)
{
  Point_set point_set;

  // Generate points on a unit sphere
  Point_generator generator(1.);
  std::size_t nb_pts = 10000;
  point_set.reserve (nb_pts);
  for (std::size_t i = 0; i < nb_pts; ++ i)
    point_set.insert (*(generator ++));


  // Add normal property and estimate normal values
  CGAL::jet_estimate_normals<CGAL::Sequential_tag> (point_set,
                                                    12); // Number of neighbors


  // Simplify point set
  CGAL::grid_simplify_point_set (point_set,
                                 0.1); // Size of grid cell

  std::vector<std::string> properties = point_set.properties();
  std::cerr << "Properties:" << std::endl;
  for (std::size_t i = 0; i < properties.size(); ++ i)
    std::cerr << " * " << properties[i] << std::endl;

  // Detect sphere with RANSAC
  Efficient_ransac ransac;
  ransac.set_input(point_set,
                   point_set.point_map(), // Call built-in property map
                   point_set.normal_map()); // Call built-in property map
  ransac.add_shape_factory<Sphere>();
  Efficient_ransac::Parameters parameters;
  parameters.probability = 0.05;
  parameters.min_points = std::size_t(point_set.size() / 3);
  parameters.epsilon = 0.01;
  parameters.cluster_epsilon = 0.5;
  parameters.normal_threshold = 0.9;   
  ransac.detect(parameters);
  
  BOOST_FOREACH(boost::shared_ptr<Efficient_ransac::Shape> shape, ransac.shapes())
    if (Sphere* sphere = dynamic_cast<Sphere*>(shape.get()))
      std::cerr << "Detected sphere of center " << sphere->center() // Center should be approx 0, 0, 0
                << " and of radius " << sphere->radius() << std::endl; // Radius should be approx 1

  return 0;
}
