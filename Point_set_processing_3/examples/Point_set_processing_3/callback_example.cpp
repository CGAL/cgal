#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Real_timer.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>

#include <boost/lexical_cast.hpp>

#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Random_points_on_sphere_3<Point> Generator;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// instance of std::function<bool(double)>
struct Progress_to_std_cerr_callback
{
  mutable std::size_t nb;
  CGAL::Real_timer timer;
  double t_start;
  mutable double t_latest;
  const std::string name;

  Progress_to_std_cerr_callback (const char* name)
    : name (name)
  {
    timer.start();
    t_start = timer.time();
    t_latest = t_start;
  }

  bool operator()(double advancement) const
  {
    // Avoid calling time() at every single iteration, which could
    // impact performances very badly
    ++ nb;
    if (advancement != 1 && nb % 100 != 0)
      return true;

    double t = timer.time();
    if (advancement == 1 || (t - t_latest) > 0.1) // Update every 1/10th of second
    {
      std::cerr << "\r" // Return at the beginning of same line and overwrite
                << name << ": " << int(advancement * 100) << "%";

      if (advancement == 1)
        std::cerr << std::endl;
      t_latest = t;
    }

    return true;
  }
};


int main (int argc, char* argv[])
{
  int N = (argc > 1) ? boost::lexical_cast<int>(argv[1]) : 1000;

  // Generate N points on a sphere of radius 100.
  std::vector<Point> points;
  points.reserve (N);
  Generator generator(100.);
  std::copy_n (generator, N, std::back_inserter(points));

  // Compute average spacing
  FT average_spacing = CGAL::compute_average_spacing<Concurrency_tag>
    (points, 6,
     CGAL::parameters::callback
     (Progress_to_std_cerr_callback("Computing average spacing")));

  // Simplify on a grid with a size of twice the average spacing
  points.erase(CGAL::grid_simplify_point_set
               (points, 2. * average_spacing,
                CGAL::parameters::callback
                (Progress_to_std_cerr_callback("Grid simplification"))),
               points.end());

  // Smooth simplified point set
  CGAL::jet_smooth_point_set<Concurrency_tag>
    (points, 6,
     CGAL::parameters::callback
     (Progress_to_std_cerr_callback("Jet smoothing")));

  return EXIT_SUCCESS;
}

