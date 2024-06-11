#include "include/test_efficient_RANSAC_generators.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

template <class K>
bool test_plane_parameters() {
  const int NB_ROUNDS = 10;
  const int NB_POINTS = 1000;

  typedef typename K::FT                                      FT;
  typedef typename CGAL::Point_with_normal_3<K>               Pwn;
  typedef typename CGAL::Vector_3<K>                          Vector;
  typedef std::vector<Pwn>                                    Pwn_vector;
  typedef typename CGAL::Identity_property_map<Pwn>           Point_map;
  typedef typename CGAL::Normal_of_point_with_normal_map<K>   Normal_map;

  typedef CGAL::Shape_detection::Efficient_RANSAC_traits<K, Pwn_vector, Point_map, Normal_map> Traits;

  typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;
  typedef CGAL::Shape_detection::Plane<Traits>            Plane;

  std::size_t success = 0;
  for (int i = 0;i<NB_ROUNDS;i++) {
    Pwn_vector points;

    FT dist = 0;
    Vector normal;
    CGAL::Bbox_3 bbox(-10, -10, -10, 10, 10, 10);

    sample_random_parallelogram_in_box(NB_POINTS, bbox, normal,
      dist, std::back_inserter(points));

    if (i >= NB_ROUNDS / 2)
      for (std::size_t j = 0;j<NB_POINTS/2;j++) {
        points.push_back(random_pwn_in<K>(bbox));
      }

    Efficient_ransac ransac;
    Traits const& traits = ransac.traits();

    ransac.template add_shape_factory<Plane>();

    ransac.set_input(points);

    // Set cluster epsilon to a high value as just the parameters of
    // the extracted primitives are to be tested.
    typename Efficient_ransac::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = 100;
    parameters.epsilon = 0.002f;
    parameters.cluster_epsilon = 1.0f;
    parameters.normal_threshold = 0.9f;

    if (!ransac.detect(parameters)) {
      std::cout << " aborted" << std::endl;
      return false;
    }

    typename Efficient_ransac::Shape_range shapes = ransac.shapes();

    if (shapes.size() != 1)
      continue;

    std::shared_ptr<Plane> pl = std::dynamic_pointer_cast<Plane>((*shapes.first));

    if (!pl)
      continue;

    const FT phi = traits.compute_scalar_product_3_object()(
      normal, pl->plane_normal());
    const FT sign = (phi < 0) ? -1.0f : 1.0f;

    const FT dist2 = pl->d();

    if (CGAL::abs(phi) < 0.98 || CGAL::abs(dist2 - sign * dist) > 0.02)
      continue;

    std::string info = pl->info();

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

  std::cout << "test_plane_parameters<CGAL::Simple_cartesian<float>> ";
  if (!test_plane_parameters<CGAL::Simple_cartesian<float> >())
    success = false;

  std::cout << "test_plane_parameters<CGAL::Simple_cartesian<double>> ";
  if (!test_plane_parameters<CGAL::Simple_cartesian<double> >())
    success = false;

  std::cout << "test_plane_parameters<CGAL::Exact_predicates_inexact_constructions_kernel> ";
  if (!test_plane_parameters<CGAL::Exact_predicates_inexact_constructions_kernel>())
    success = false;

  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}
