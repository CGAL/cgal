#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/read_points.h>
#include <CGAL/estimate_scale.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/grid_simplify_point_set.h>

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <vector>
#include <fstream>

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;

int main (int argc, char** argv)
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("points_3/sphere_20k.xyz");

  CGAL::Timer task_timer;

  // read input
  std::vector<Point_3> points;
  if(!CGAL::IO::read_points(fname, std::back_inserter(points)))
  {
    std::cerr << "Error: can't read input file" << std::endl;
    return EXIT_FAILURE;
  }

  // estimate k scale
  task_timer.start();
  std::size_t k_scale = CGAL::estimate_global_k_neighbor_scale (points);
  task_timer.stop();

  // Example: use estimated k as scale for jet smoothing
  CGAL::jet_smooth_point_set<Concurrency_tag>(points, static_cast<unsigned int>(k_scale));

  // estimate range scale
  task_timer.start();
  FT range_scale = CGAL::estimate_global_range_scale (points);
  task_timer.stop();

  // Example: use estimated range for grid simplification
  points.erase(CGAL::grid_simplify_point_set(points, range_scale), points.end());

  // print some informations on runtime
  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  double time = task_timer.time();

  std::cout << "Scales computed in " << time << " second(s) using "
            << (memory>>20) << " MiB of memory:" << std::endl;
  std::cout << " * Global K scale: " << k_scale << std::endl;
  std::cout << " * Global range scale: " << range_scale << std::endl;

  return EXIT_SUCCESS;
}

