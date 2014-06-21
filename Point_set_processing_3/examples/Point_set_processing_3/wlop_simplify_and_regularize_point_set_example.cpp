#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
  const std::string INPUT_FILENAME_WITHOUT_EXT = "data/saint_jean_370K";
  const std::string INPUT_FILENAME = INPUT_FILENAME_WITHOUT_EXT + ".xyz";
  const std::string OUTPUT_FILENAME = INPUT_FILENAME_WITHOUT_EXT + "_WLOPED.xyz";

  // Reads a .xyz point set file in points[], *with normals*.
	// TOFIX: do we have an extension for points with normals?
  std::vector<Point> points;
  std::ifstream stream(INPUT_FILENAME.c_str());

  if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file "
              << INPUT_FILENAME.c_str()  << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Point> output;

  //with default parameters begins
  CGAL::wlop_simplify_and_regularize_point_set
                          <CGAL::Parallel_tag> // parallel version
                          (points.begin(), 
                           points.end(),
                           std::back_inserter(output));
  //with default parameters ends

  //with all parameters begins
  /*
  //Algorithm parameters
  const double select_percentage = 2;       // percentage of points to retain.
  const double neighbor_radius = 0.03;      // neighbors size.
  const unsigned int iter_number = 30;      // number of iterations.
  const bool require_uniform_sampling = true;  // Optional pre-processing. 
                                       
  CGAL::wlop_simplify_and_regularize_point_set
                          <CGAL::Parallel_tag> // parallel version
                          (points.begin(),
                           points.end(),
                           back_inserter(output),
                           select_percentage, 
                           neighbor_radius,
                           iter_number,
                           require_uniform_sampling);
  */
  //with all parameters ends
  
  std::ofstream out(OUTPUT_FILENAME.c_str()); 
  if (!out || !CGAL::write_xyz_points(
      out, output.begin(), output.end()))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}




