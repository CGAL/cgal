#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_off_points.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Kernel> Point_set;

int main (int argc, char** argv)
{
  std::ifstream f (argc > 1 ? argv[1] : "data/oni.xyz");

  Point_set point_set;
  point_set.add_normal_property();

  // Reading input in XYZ format
  if (!f ||
      !CGAL::read_xyz_points_and_normals
      (f,
       point_set.index_back_inserter(),
       point_set.point_push_pmap(), // Use push property map for creating new items
       point_set.normal_push_pmap())) // Same
    {
      std::cerr << "Can't read input file " << std::endl;
      return EXIT_FAILURE;
    }

  // Normalization + inversion of normal vectors
  for (Point_set::iterator it = point_set.begin(); it != point_set.end(); ++ it)
    {
      Vector n = point_set.normal(it);
      n = - n / std::sqrt (n * n);
      point_set.normal(it) = n;
    }

  // Writing result in OFF format
  std::ofstream out("normalized_normals.off");
  if (!out ||
      !CGAL::write_off_points_and_normals
      (out, point_set.begin(), point_set.end(),
       point_set.point_pmap(), point_set.normal_pmap())) // Use regular property map for accessing
    {
      return EXIT_FAILURE;
    }
    
  return EXIT_SUCCESS;
}
