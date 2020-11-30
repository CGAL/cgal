#include <CGAL/Simple_cartesian.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>

#include <CGAL/IO/write_ply_points.h>

#include <CGAL/Real_timer.h>

#include <boost/function_output_iterator.hpp>

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

void test_random_planes(std::size_t nb_planes);

int main()
{
  test_random_planes(1);
  test_random_planes(10);
  test_random_planes(100);
  test_random_planes(1000);
  test_random_planes(10000);

  return EXIT_SUCCESS;
}

void test_random_planes(std::size_t nb_planes)
{
  CGAL::Random random;
  std::cerr << "[TEST ON " << nb_planes << " RANDOM PLANES, SEED = "
            << random.get_seed() << "] ";

  std::vector<double> time_rg;
  std::vector<std::size_t> nb_detected_rg;
  std::vector<double> nb_unassigned_rg;
  std::vector<double> time_ransac;
  std::vector<std::size_t> nb_detected_ransac;
  std::vector<double> nb_unassigned_ransac;

  CGAL::Real_timer timer;
  double timeout = 120.; // 2 minutes timeout

  timer.start();
  std::size_t nb_runs = 500;
  for (std::size_t run = 0; run < nb_runs; ++ run)
  {
    std::size_t min_points = random.get_int(10, 200);
    std::size_t max_points = random.get_int(200, 10000);
    double cluster_epsilon = random.get_double(0.01, 10.);
    double epsilon = random.get_double (0., cluster_epsilon / 10.);
    double normal_threshold = random.get_double (0.75, 0.99);

    double domain_size = cluster_epsilon * std::sqrt(nb_planes) * min_points;
    double spacing = 0.8 * cluster_epsilon / std::sqrt(2); // smaller than diagonal
    double noise = 0.5 * epsilon;

    Point_set points;
    for (std::size_t i = 0; i < nb_planes; ++ i)
    {
      // Generate random plane
      Point_3 origin (random.get_double(-domain_size, domain_size),
                      random.get_double(-domain_size, domain_size),
                      random.get_double(-domain_size, domain_size));
      Vector_3 base1 (random.get_double(-1,1), random.get_double(-1,1), random.get_double(-1,1));
      base1 = base1 / std::sqrt(base1 * base1);
      Vector_3 base2 (random.get_double(-1,1), random.get_double(-1,1), random.get_double(-1,1));
      base2 = base2 / std::sqrt(base2 * base2);

      Vector_3 normal = CGAL::cross_product (base1, base2);
      normal = normal / std::sqrt (normal * normal);

      std::size_t nb_points = random.get_int (min_points, max_points);

      std::size_t nb_x = std::size_t(std::sqrt (double(nb_points)));
      std::size_t nb_y = std::size_t(nb_points / double(nb_x)) + 1;

      for (std::size_t j = 0; j < nb_x; ++ j)
        for (std::size_t k = 0; k < nb_y; ++ k)
          points.emplace_back (origin + j * spacing * base1 + k * spacing * base2
                               + normal * random.get_double(-noise/2, noise/2),
                               normal);
    }

    // std::ofstream test("dump.pwn");
    // for (const auto& p : points)
    //   test << p.first << " " << p.second << std::endl;
    CGAL::Real_timer t;
    t.start();
    RG_query rg_query (points, cluster_epsilon);
    RG_region rg_region (points, epsilon, normal_threshold, min_points);
    Region_growing region_growing (points, rg_query, rg_region);
    std::size_t nb_detected = 0;
    std::size_t nb_unassigned = 0;
    region_growing.detect (boost::make_function_output_iterator ([&](const auto&) { ++ nb_detected; }));
    region_growing.unassigned_items (boost::make_function_output_iterator ([&](const auto&) { ++ nb_unassigned; }));
    t.stop();

    time_rg.push_back (t.time());
    nb_detected_rg.push_back (nb_detected);
    nb_unassigned_rg.push_back (nb_unassigned / double(points.size()));
    t.reset();

    t.start();
    RANSAC ransac;
    ransac.template add_shape_factory<RANSAC_plane>();
    ransac.set_input(points);
    typename RANSAC::Parameters parameters;
    parameters.probability = 0.05f;
    parameters.min_points = min_points;
    parameters.epsilon = epsilon;
    parameters.cluster_epsilon = cluster_epsilon;
    parameters.normal_threshold = normal_threshold;
    ransac.detect(parameters);
    t.stop();

    time_ransac.push_back (t.time());
    nb_detected_ransac.push_back (ransac.shapes().size());
    nb_unassigned_ransac.push_back (ransac.number_of_unassigned_points() / double(points.size()));

#if 0
    if (ransac.shapes().size() == 0)
    {
      std::cerr << "Detected 0 shapes with " << std::endl
                << " * min points = " << min_points << std::endl
                << " * epsilon = " << epsilon << std::endl
                << " * cluster_epsilon = " << cluster_epsilon << std::endl
                << " * normal_threshold = " << normal_threshold << std::endl;

      std::ofstream ofile("0shapes.ply", std::ios::binary);
//      CGAL::set_binary_mode (ofile);
      CGAL::write_ply_points (ofile, points,
                              CGAL::parameters::point_map(Point_map()).
                              normal_map(Normal_map()));

      exit(0);
    }
#endif

    if (timer.time() > timeout)
    {
      nb_runs = run + 1;
      break;
    }
  }

  std::cerr << "on " << nb_runs << " runs" << std::endl;


  std::sort (time_ransac.begin(), time_ransac.end());
  std::sort (nb_detected_ransac.begin(), nb_detected_ransac.end());
  std::sort (nb_unassigned_ransac.begin(), nb_unassigned_ransac.end());
  std::sort (time_rg.begin(), time_rg.end());
  std::sort (nb_detected_rg.begin(), nb_detected_rg.end());
  std::sort (nb_unassigned_rg.begin(), nb_unassigned_rg.end());

  std::cerr << " * Region Growing" << std::endl
            << "   - took between " << time_rg.front() << "s and " << time_rg.back()
            << "s, median = " << time_rg[nb_runs / 2] << "s" << std::endl
            << "   - detected between " << nb_detected_rg.front() << " and "
            << nb_detected_rg.back() << " planes (" << 100. * nb_detected_rg.front() / double(nb_planes)
            << "% to " << 100. * nb_detected_rg.back() / double(nb_planes) << "%), median = "
            << nb_detected_rg[nb_runs / 2] << " planes (" << 100. * nb_detected_rg[nb_runs / 2] / double(nb_planes)
            << "%)" << std::endl
            << "   - left between " << 100. * nb_unassigned_rg.front()
            << "% and " << 100. * nb_unassigned_rg.back() << "% of unassigned points, median = "
            << 100. * nb_unassigned_rg[nb_runs / 2] << "%" << std::endl;
  std::cerr << " * Efficient RANSAC" << std::endl
            << "   - took between " << time_ransac.front() << "s and " << time_ransac.back()
            << "s, median = " << time_ransac[nb_runs / 2] << "s" << std::endl
            << "   - detected between " << nb_detected_ransac.front() << " and "
            << nb_detected_ransac.back() << " planes (" << 100. * nb_detected_ransac.front() / double(nb_planes)
            << "% to " << 100. * nb_detected_ransac.back() / double(nb_planes) << "%), median = "
            << nb_detected_ransac[nb_runs / 2] << " planes (" << 100. * nb_detected_ransac[nb_runs / 2] / double(nb_planes)
            << "%)" << std::endl
            << "   - left between " << 100. * nb_unassigned_ransac.front()
            << "% and " << 100. * nb_unassigned_ransac.back() << "% of unassigned points, median = "
            << 100. * nb_unassigned_ransac[nb_runs / 2] << "%" << std::endl;
}
