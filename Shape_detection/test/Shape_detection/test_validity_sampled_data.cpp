#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

#ifdef CGAL_TEST_RANSAC_PROTOTYPE
#include <RansacShapeDetector.h>
#include <CylinderPrimitiveShapeConstructor.h>
#include <PlanePrimitiveShapeConstructor.h>
#include <CylinderPrimitiveShape.h>
#endif

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_ply_points.h>

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
using Point_map = CGAL::First_of_pair_property_map<Pwn>;
using Normal_map = CGAL::Second_of_pair_property_map<Pwn>;

using RG_query = SD::Point_set::Sphere_neighbor_query<Kernel, Point_set, Point_map>;
using RG_region = SD::Point_set::Least_squares_plane_fit_region<Kernel, Point_set, Point_map, Normal_map>;
using Region_growing = SD::Region_growing<Point_set, RG_query, RG_region>;

using RANSAC_traits = SD::Efficient_RANSAC_traits<Kernel, Point_set, Point_map, Normal_map>;
using RANSAC = SD::Efficient_RANSAC<RANSAC_traits>;
using RANSAC_plane = SD::Plane<RANSAC_traits>;

void test_copied_point_cloud (const Point_set& points, std::size_t nb);

int main (int argc, char** argv)
{
  Point_set points;
  const char* ifilename = (argc > 1) ? argv[1] : "data/point_set_3.xyz";
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
     CGAL::make_transform_iterator_from_property_map (points.end(), Point_map()));

  typename RANSAC::Parameters parameters;
  parameters.probability = 0.01;
  parameters.min_points = 100;
  parameters.epsilon = 0.01;
  parameters.cluster_epsilon = 0.01;
  parameters.normal_threshold = 0.95;

  CGAL::Real_timer t;
  t.start();
  RG_query rg_query (points, parameters.cluster_epsilon);
  RG_region rg_region (points, parameters.epsilon, parameters.normal_threshold, parameters.min_points);
  Region_growing region_growing (points, rg_query, rg_region);
  std::size_t nb_detected = 0;
  std::size_t nb_unassigned = 0;
  region_growing.detect (boost::make_function_output_iterator ([&](const auto&) { ++ nb_detected; }));
  region_growing.unassigned_items (boost::make_function_output_iterator ([&](const auto&) { ++ nb_unassigned; }));
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

  // RANSAC should at least detect 75% of shapes
  assert (detected_ransac[detected_ransac.size() / 2] > std::size_t(0.75 * ground_truth));

#ifdef CGAL_TEST_RANSAC_PROTOTYPE
  {
    CGAL::Real_timer timer;
    double timeout = 120.; // 2 minute timeout
    timer.start();
    std::size_t nb_runs = 500;
    std::vector<std::size_t> detected_ransac;
    std::vector<double> times_ransac;
    for (std::size_t run = 0; run < nb_runs; ++ run)
    {
      PointCloud proto_points;
      proto_points.reserve (points.size());

      Point Pt;
      for (std::size_t i = 0; i < points.size(); ++i)
      {
        Pt.pos[0] = static_cast<float>(points[i].first.x());
        Pt.pos[1] = static_cast<float>(points[i].first.y());
        Pt.pos[2] = static_cast<float>(points[i].first.z());
        Pt.normal[0] = static_cast<float>(points[i].second.x());
        Pt.normal[1] = static_cast<float>(points[i].second.y());
        Pt.normal[2] = static_cast<float>(points[i].second.z());
#ifdef POINTSWITHINDEX
        Pt.index = i;
#endif
        proto_points.push_back(Pt);
      }

      //manually set bounding box!
      Vec3f cbbMin, cbbMax;
      cbbMin[0] = static_cast<float>(bbox.xmin());
      cbbMin[1] = static_cast<float>(bbox.ymin());
      cbbMin[2] = static_cast<float>(bbox.zmin());
      cbbMax[0] = static_cast<float>(bbox.xmax());
      cbbMax[1] = static_cast<float>(bbox.ymax());
      cbbMax[2] = static_cast<float>(bbox.zmax());
      proto_points.setBBox(cbbMin, cbbMax);

      // Sets parameters for shape detection.
      RansacShapeDetector::Options options;
      options.m_epsilon = parameters.epsilon;
      options.m_bitmapEpsilon = parameters.cluster_epsilon;
      options.m_normalThresh = parameters.normal_threshold;
      options.m_probability = parameters.probability;
      options.m_minSupport = parameters.min_points;

      CGAL::Real_timer t;
      t.start();
      RansacShapeDetector ransac (options); // the detector object
      ransac.Add (new PlanePrimitiveShapeConstructor());
      MiscLib::Vector<std::pair<MiscLib::RefCountPtr<PrimitiveShape>, std::size_t> > shapes; // stores the detected shapes
      ransac.Detect (proto_points, 0, proto_points.size(), &shapes);
      t.stop();

      detected_ransac.emplace_back (shapes.size());
      times_ransac.emplace_back (t.time() * 1000);
      if (timer.time() > timeout)
      {
        nb_runs = run + 1;
        break;
      }
    }

    std::sort (detected_ransac.begin(), detected_ransac.end());
    std::sort (times_ransac.begin(), times_ransac.end());
    std::cerr << "RANSAC (proto) = " << detected_ransac[detected_ransac.size() / 2]
              << " planes (" << times_ransac[times_ransac.size() / 2] << "ms)     (on "
              << nb_runs << " runs, planes["
              << detected_ransac.front() << ";" << detected_ransac.back() << "], time["
              << times_ransac.front() << ";" << times_ransac.back() << "])" << std::endl;
  }
#endif

  std::cerr << std::endl;
}
