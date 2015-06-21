#include <CGAL/Simple_cartesian.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

int main(int argc, char** argv)
{
  const char* input_filename = (argc>1)?argv[1]:"data/sphere_20k.xyz";
  const char* output_filename = (argc>2)?argv[2]:"data/sphere_20k_WLOPED.xyz";

  // Reads a .xyz point set file in points[]
  std::vector<Point> points;
  std::ifstream stream(input_filename);

  if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << input_filename  << std::endl;

    return EXIT_FAILURE;
  }

  std::vector<Point> output;

  //parameters
  const double retain_percentage = 2;   // percentage of points to retain.
  const double neighbor_radius = 0.5;   // neighbors size.

  CGAL::wlop_simplify_and_regularize_point_set
                          <CGAL::Parallel_tag> // parallel version
                          (points.begin(), 
                           points.end(),
                           std::back_inserter(output),
                           retain_percentage,
                           neighbor_radius
                           );
  
  std::ofstream out(output_filename);
  if (!out || !CGAL::write_xyz_points(
      out, output.begin(), output.end()))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}




