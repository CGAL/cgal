#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/regularize_and_simplify_point_set.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/tags.h>

#include <vector>
#include <fstream>
#include <time.h>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
  const std::string INPUT_FILENAME_WITHOUT_EXT = "data/sphere_20k";
  //const std::string INPUT_FILENAME_WITHOUT_EXT = "data/saint_jean_370K";

  // Reads a .xyz point set file in points[], *with normals*.
  std::vector<Point> points;
  std::ifstream stream(INPUT_FILENAME_WITHOUT_EXT + ".xyz");

  if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file "
              << INPUT_FILENAME_WITHOUT_EXT << ".xyz" << std::endl;
    return EXIT_FAILURE;
  }

  //Algorithm parameters
  const double retain_percentage = 5;   // percentage of points to retain.
  const double neighbor_radius = 0.25;   // neighbors size.
  const unsigned int iter_number = 35;     // number of iterations.
  const bool need_compute_density = true;  // if needed to compute density.

  // Make room for sample points
  std::vector<Point> points_sampled;
  points_sampled.resize(points.size() * (retain_percentage / 100.));

  CGAL::Timer task_timer;
  task_timer.start();
  std::cout << "Run algorithm example: " << std::endl;

  // Run algorithm 
  //std::vector<Point>::const_iterator sample_points_begin =
  //  CGAL::wlop_simplify_and_regularize_point_set(
  //          points.begin(), 
  //          points.end(), 
  //          retain_percentage, 
  //          neighbor_radius,
  //          iter_number,
  //          need_compute_density);

  std::vector<Point>::const_iterator sample_points_begin =
  CGAL::wlop_simplify_and_regularize_point_set<CGAL::Parallel_tag>(
            points.begin(), 
            points.end());

  // Copy results to sample points
  std::copy(sample_points_begin,
            static_cast<std::vector<Point>::const_iterator>(points.end()),
            points_sampled.begin());


  long memory = CGAL::Memory_sizer().virtual_size();
  std::cout << "total done: " << task_timer.time() << " seconds, " 
            << (memory>>20) << " Mb allocated" << std::endl;

  // Save point set.
  std::ofstream out(INPUT_FILENAME_WITHOUT_EXT + "_WLOPED.xyz");  
  if (!out ||
     !CGAL::write_xyz_points(
      out, points_sampled.begin(), points_sampled.end()))
  {
    return EXIT_FAILURE;
  }

  system("Pause");

  return EXIT_SUCCESS;
}




