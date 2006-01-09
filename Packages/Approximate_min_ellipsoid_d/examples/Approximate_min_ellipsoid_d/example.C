#include <vector>
#include <iostream>

#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>

typedef CGAL::Cartesian_d<double>                              Kernel;
typedef CGAL::MP_Float                                         ET;
typedef CGAL::Approximate_min_ellipsoid_d_traits_d<Kernel, ET> Traits;
typedef Traits::Point                                          Point;
typedef std::vector<Point>                                     Point_list;
typedef CGAL::Approximate_min_ellipsoid_d<Traits>              AME;

int main()
{
  const int      n = 1000;                // number of points
  const int      d = 2;                   // dimension
  const double eps = 0.01;                // approximation ratio is (1+eps)

  // create a set of random points:
  Point_list P;
  CGAL::Random_points_in_iso_box_d<Point> rpg(d,100.0);
  for (int i = 0; i < n; ++i) {
    P.push_back(*rpg);
    ++rpg;
  }

  // compute approximation:
  Traits traits;
  AME mel(eps, P.begin(), P.end(), traits);
  
  // write EPS file:
  if (mel.is_full_dimensional() && d == 2)
    mel.write_eps("example.eps");
}
