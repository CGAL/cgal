#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/grid_simplify_point_set.h>

#include <CGAL/Shape_detection_3.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;

typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits
<Kernel, Point_set, Point_set::Point_pmap, Point_set::Vector_pmap> Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>       Efficient_ransac;
typedef CGAL::Shape_detection_3::Sphere<Traits>                 Sphere;

int main (int, char**)
{
  Point_set point_set;

  // Generate points on a sphere
  std::size_t nb_pts = 10000;
  point_set.reserve (nb_pts);
  for (std::size_t i = 0; i < nb_pts; ++ i)
    {
      double sintheta = 2 * rand () / (double)RAND_MAX;
      double theta = std::asin (sintheta);
      double phi = 2 * M_PI * (rand () / (double)RAND_MAX) - M_PI;
      Point p (std::cos (theta) * std::cos (phi),
      	       std::cos (theta) * std::sin (phi),
	       std::sin (theta));
      point_set.push_back (p);
    }

  // Add normal property and estimate normal values
  point_set.add_normal_property();
  CGAL::jet_estimate_normals<CGAL::Sequential_tag> (point_set.begin(), point_set.end(),
                                                    point_set.point_pmap(),
                                                    point_set.normal_pmap(),
                                                    12); // Number of neighbors


  // Simplify point set
  point_set.remove_from
    (CGAL::grid_simplify_point_set (point_set.begin(), point_set.end(),
                                    point_set.point_pmap(), 
                                    0.1)); // Size of grid cell

  std::cerr << point_set.properties();

  // Detect sphere with RANSAC
  Efficient_ransac ransac;
  ransac.set_input(point_set, point_set.point_pmap(), point_set.normal_pmap());
  ransac.add_shape_factory<Sphere>();
  Efficient_ransac::Parameters parameters;
  parameters.probability = 0.05;
  parameters.min_points = point_set.size() / 3.;
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
