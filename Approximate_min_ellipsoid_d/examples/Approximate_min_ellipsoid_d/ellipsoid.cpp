#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>

#include <vector>
#include <iostream>

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
  CGAL::Random_points_in_cube_d<Point> rpg(d,100.0);
  for (int i = 0; i < n; ++i) {
    P.push_back(*rpg);
    ++rpg;
  }

  // compute approximation:
  Traits traits;
  AME ame(eps, P.begin(), P.end(), traits);

  // write EPS file:
  if (ame.is_full_dimensional() && d == 2)
    ame.write_eps("example.eps");

  // output center coordinates:
  std::cout << "Cartesian center coordinates: ";
  for (AME::Center_coordinate_iterator c_it = ame.center_cartesian_begin();
       c_it != ame.center_cartesian_end();
       ++c_it)
    std::cout << *c_it << ' ';
  std::cout << ".\n";

  if (d == 2 || d == 3) {
    // output  axes:
    AME::Axes_lengths_iterator axes = ame.axes_lengths_begin();
    for (int i = 0; i < d; ++i) {
      std::cout << "Semiaxis " << i << " has length " << *axes++  << "\n"
                << "and Cartesian coordinates ";
      for (AME::Axes_direction_coordinate_iterator
             d_it = ame.axis_direction_cartesian_begin(i);
           d_it != ame.axis_direction_cartesian_end(i); ++d_it)
        std::cout << *d_it << ' ';
      std::cout << ".\n";
    }
  }
}
