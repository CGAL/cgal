#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/jet_estimate_normals.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Kernel> Point_set;

int main (int, char**)
{
  Point_set point_set;

  // Generate points on a sphere
  std::size_t nb_pts = 10000;
  point_set.reserve (nb_pts);
  for (std::size_t i = 0; i < nb_pts; ++ i)
    {
      double sintheta = 2 * rand () / (double)RAND_MAX;
      double theta = std::asin (sintheta);
      double phi = 2 * M_PI * (rand () / (double)RAND_MAX) - M_PI;
      Point p (std::cos (theta) * std::cos (phi),
      	       std::cos (theta) * std::sin (phi),
	       std::sin (theta));
      point_set.push_back (p);
    }

  point_set.add_normal_property();


  CGAL::jet_estimate_normals<CGAL::Sequential_tag> (point_set.begin(), point_set.end(),
                                                    point_set.point_pmap(),
                                                    point_set.normal_pmap(),
                                                    12); // Number of neighbors
  

  return 0;
}
