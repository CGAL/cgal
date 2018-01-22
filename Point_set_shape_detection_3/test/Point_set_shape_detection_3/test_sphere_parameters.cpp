#include "generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>


template <class K>
bool test_sphere_parameters() {
  const int NB_ROUNDS = 10;
  const int NB_POINTS = 1000;
  
  typedef typename K::FT                                      FT;
  typedef CGAL::Point_with_normal_3<K>                        Pwn;
  typedef CGAL::Point_3<K>                                    Point;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef CGAL::Identity_property_map<Pwn>                    Point_map;
  typedef CGAL::Normal_of_point_with_normal_pmap<K>           Normal_map;

  typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<
    K, Pwn_vector, Point_map, Normal_map>                     Traits;

  typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>   Efficient_ransac;
  typedef CGAL::Shape_detection_3::Sphere<Traits>             Sphere;

  std::size_t success = 0;

  for (int i = 0;i<NB_ROUNDS;i++) {
    Pwn_vector points;

    // generate random points on random sphere
    FT radius = 0;
    Point center;
    CGAL::Bbox_3 bbox(-10, -10, -10, 10, 10, 10);

    sample_random_sphere_in_box(NB_POINTS, bbox, center,
      radius, std::back_inserter(points));

    // add outliers in second half of rounds
    if (i >= NB_ROUNDS / 2)
      for (std::size_t j = 0; j < NB_POINTS / 2; j++) 
        points.push_back(random_pwn_in<K>(bbox));

    Efficient_ransac ransac;
    Traits traits = ransac.traits();

    ransac.template add_shape_factory<Sphere>();

    ransac.set_input(points);

    // Set cluster epsilon to a high value as just the parameters of
    // the extracted primitives are to be tested.
    typename Efficient_ransac::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = NB_POINTS/10;
    parameters.epsilon = 0.002f;
    parameters.cluster_epsilon = 1.0f;
    parameters.normal_threshold = 0.9f;

    if (!ransac.detect(parameters)) {
      std::cout << " aborted" << std::endl;
      return false;
    }

    typename Efficient_ransac::Shape_range shapes = ransac.shapes();

    // check: unique shape detected
    if (shapes.size() != 1)
      continue;

    boost::shared_ptr<Sphere> sphere = boost::dynamic_pointer_cast<Sphere>((*shapes.first));

    // check: shape detected is a cylinder
    if (!sphere)
      continue;

    // Check radius and alignment with axis.
    if (CGAL::abs(radius - sphere->radius()) > (FT) 0.02)
      continue;

    // Check center.
    Point pos = sphere->center();
    FT center_pos_sqlen = traits.compute_squared_length_3_object()(
      traits.construct_vector_3_object()(center, pos));
    if (center_pos_sqlen > FT(0.0004))
      continue;

    std::string info = sphere->info();

    success++;
  }

  if (success >= NB_ROUNDS * 0.8) {
    std::cout << " succeeded" << std::endl;
    return true;
  }
  else {
    std::cout << " failed" << std::endl;
    return false;
  }
}


int main() {
  bool success = true;

  std::cout << "test_sphere_parameters<CGAL::Simple_cartesian<float>> ";
  if (!test_sphere_parameters<CGAL::Simple_cartesian<float> >()) 
    success = false;

  std::cout << "test_sphere_parameters<CGAL::Simple_cartesian<double>> ";
  if (!test_sphere_parameters<CGAL::Simple_cartesian<double> >())
    success = false;

  std::cout << "test_sphere_parameters<CGAL::Exact_predicates_inexact_constructions_kernel> ";
  if (!test_sphere_parameters<CGAL::Exact_predicates_inexact_constructions_kernel>()) 
    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
