#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(int argc, char*argv[])
{
  // Reads a .xyz point set file in points[].
  std::vector<Point> points;
  const char* fname = (argc>1)?argv[1]:"data/oni.xyz";
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Read " << points.size () << " point(s)" << std::endl;

  CGAL::Timer task_timer; task_timer.start();

  // simplification by clustering using erase-remove idiom
  points.erase (CGAL::hierarchy_simplify_point_set (points,
                                                    CGAL::parameters::size(100). // Max cluster size
                                                    maximum_variation(0.01)), // Max surface variation
                points.end ());

  std::size_t memory = CGAL::Memory_sizer().virtual_size();

  std::cout << points.size () << " point(s) kept, computed in "
            << task_timer.time() << " seconds, "
            << (memory>>20) << " Mib allocated." << std::endl;

  std::ofstream f ("out.xyz");
  f.precision(17);
  CGAL::write_xyz_points (f, points);

  return EXIT_SUCCESS;
}

