#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>


#include <CGAL/IO/XYZ.h>
#include <CGAL/IO/PLY.h>

#include <CGAL/Real_timer.h>

#include <boost/iterator/function_output_iterator.hpp>

// Uncomment this line to run expensive test
// #define CGAL_SHAPE_DETECTION_RUN_EXPENSIVE_TESTS

namespace SD = CGAL::Shape_detection;

using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;

using Pwn = std::pair<Point_3, Vector_3>;
using Point_set = std::vector<Pwn>;
using Item           = Point_set::const_iterator;

using Deref_map         = CGAL::Dereference_property_map<const Pwn, Item>;
using Point_map = CGAL::First_of_pair_property_map<Pwn>;
using Normal_map = CGAL::Second_of_pair_property_map<Pwn>;
using RG_Point_map         = CGAL::Compose_property_map<Deref_map, Point_map>;
using RG_Normal_map        = CGAL::Compose_property_map<Deref_map, Normal_map>;

using RG_query = SD::Point_set::Sphere_neighbor_query<Kernel, Item, RG_Point_map>;
using RG_region = SD::Point_set::Least_squares_plane_fit_region<Kernel, Item, RG_Point_map, RG_Normal_map>;
using Region_growing = SD::Region_growing<RG_query, RG_region>;

using RANSAC_traits = SD::Efficient_RANSAC_traits<Kernel, Point_set, Point_map, Normal_map>;
using RANSAC = SD::Efficient_RANSAC<RANSAC_traits>;
using RANSAC_plane = SD::Plane<RANSAC_traits>;

void test_copied_point_cloud (const Point_set& points, std::size_t nb);

int main (int argc, char** argv)
{
  Point_set points;
  const std::string ifilename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/building.xyz");
  std::ifstream ifile(ifilename);

  if (!ifile ||
      !CGAL::IO::read_XYZ(
      ifile,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map())))
  {
    std::cerr << "Reading error" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::cpp98::random_shuffle (points.begin(), points.end());

  test_copied_point_cloud (points, 1);
  test_copied_point_cloud (points, 2);
#ifdef CGAL_SHAPE_DETECTION_RUN_EXPENSIVE_TESTS
  test_copied_point_cloud (points, 5);
  test_copied_point_cloud (points, 10);
  test_copied_point_cloud (points, 20);
  test_copied_point_cloud (points, 50);
#endif

  return EXIT_SUCCESS;
}

void test_copied_point_cloud (const Point_set& original_points, std::size_t nb)
{
  CGAL::Bbox_3 bbox = CGAL::bbox_3
    (CGAL::make_transform_iterator_from_property_map (original_points.begin(), Point_map()),
     CGAL::make_transform_iterator_from_property_map (original_points.end(), Point_map()));

  std::size_t ground_truth = 6*nb*nb+1;
  std::cerr << "Ground truth   = " << ground_truth << " planes" << std::endl;

  Point_set points;
  points.reserve (nb * nb * original_points.size());
  for (std::size_t x = 0; x < nb; ++ x)
    for (std::size_t y = 0; y < nb; ++ y)
    {
      Vector_3 shift ((bbox.xmax() - bbox.xmin()) * x, (bbox.ymax() - bbox.ymin()) * y, 0);
      for (std::size_t i = 0; i < original_points.size(); ++ i)
      {
        Vector_3 noise (CGAL::get_default_random().get_double(-0.0001, 0.0001),
                        CGAL::get_default_random().get_double(-0.0001, 0.0001),
                        CGAL::get_default_random().get_double(-0.0001, 0.0001));
        points.emplace_back (original_points[i].first + shift + noise,
                             original_points[i].second);
      }
    }

  bbox = CGAL::bbox_3
    (CGAL::make_transform_iterator_from_property_map (points.begin(), Point_map()),
    CGAL::make_transform_iterator_from_property_map(points.end(), Point_map()));

  typename RANSAC::Parameters parameters;
  parameters.probability = 0.01;
  parameters.min_points = 100;
  parameters.epsilon = 0.01;
  parameters.cluster_epsilon = 0.01;
  parameters.normal_threshold = 0.95;

  CGAL::Real_timer t;
  t.start();
  RG_query rg_query (
    points,
    CGAL::parameters::sphere_radius(parameters.cluster_epsilon));
  RG_region rg_region (
    CGAL::parameters::
    maximum_distance(parameters.epsilon).
    cosine_of_maximum_angle(parameters.normal_threshold).
    minimum_region_size(parameters.min_points));
  Region_growing region_growing (points, rg_query, rg_region);
  std::size_t nb_detected = 0;
  std::size_t nb_unassigned = 0;
  region_growing.detect (boost::make_function_output_iterator ([&](const auto&) { ++ nb_detected; }));
  region_growing.unassigned_items(points, boost::make_function_output_iterator ([&](const auto&) { ++ nb_unassigned; }));
  t.stop();
  std::cerr << "Region Growing = " << nb_detected << " planes (" << 1000 * t.time() << "ms)" << std::endl;

  assert (nb_detected == ground_truth);

#ifdef CGAL_SHAPE_DETECTION_RUN_EXPENSIVE_TESTS
  double timeout = 120; // 2 minutes timeout
  std::size_t nb_runs = 500;
#else
  double timeout = 60; // 1 minute timeout
  std::size_t nb_runs = 20; //
#endif

  CGAL::Real_timer timer;
  timer.start();

  std::vector<std::size_t> detected_ransac;
  std::vector<double> times_ransac;
  for (std::size_t run = 0; run < nb_runs; ++ run)
  {
    Point_set ordered_points = points;

    CGAL::Real_timer t;
    t.start();
    RANSAC ransac;
    ransac.template add_shape_factory<RANSAC_plane>();
    ransac.set_input(ordered_points);
    ransac.detect(parameters);
    t.stop();

    detected_ransac.emplace_back (ransac.shapes().size());
    times_ransac.emplace_back (t.time() * 1000);
    if (timer.time() > timeout)
    {
      nb_runs = run + 1;
      break;
    }
  }

  std::sort (detected_ransac.begin(), detected_ransac.end());
  std::sort (times_ransac.begin(), times_ransac.end());
  std::cerr << "RANSAC         = " << detected_ransac[detected_ransac.size() / 2]
            << " planes (" << times_ransac[times_ransac.size() / 2] << "ms)     (on "
            << nb_runs << " runs, planes["
            << detected_ransac.front() << ";" << detected_ransac.back() << "], time["
            << times_ransac.front() << ";" << times_ransac.back() << "])" << std::endl;

  // RANSAC should detect at least 75% of shapes.
  assert (detected_ransac[detected_ransac.size() / 2] > std::size_t(0.75 * ground_truth));

  std::cerr << std::endl;
}
