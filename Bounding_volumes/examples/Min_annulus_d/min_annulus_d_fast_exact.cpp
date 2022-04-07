// computes the smallest enclosing annulus of two point
// sets on nested squares in R^2,  using double
// as input type and some internal EXACT floating point type;
// the fast type double is also safely used for many of the
// internal computations
#include <CGAL/Min_annulus_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_2.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_integer.h>
typedef CGAL::Exact_integer ET;

#include <iostream>
#include <cassert>

// use an inexact kernel...
typedef CGAL::Homogeneous<double>                          K;
typedef K::Point_2                                         Point;
// ... and the EXACT traits class based on the inexcat kernel
typedef CGAL::Min_sphere_annulus_d_traits_2<K, ET, double> Traits;
typedef CGAL::Min_annulus_d<Traits>                        Min_annulus;



int main()
{
  // points on the squares [-1,1]^2 and [-2,2]^2
  Point P[8] = { Point(-1,-1), Point(-1,1), Point(1,-1), Point(1,1),
                 Point(-2,-2), Point(-2,2), Point(2,-2), Point(2,2)};

  Min_annulus ma(P, P+8);
  assert (ma.is_valid());

  // get center of annulus
  Min_annulus::Coordinate_iterator coord_it;

  std::cout << "center:"; // homogeneous point, (0,0,1)
  for (coord_it = ma.center_coordinates_begin();
       coord_it != ma.center_coordinates_end();
       ++coord_it)
    std::cout << " " << CGAL::to_double(*coord_it);
  std::cout << std::endl;

  // get inner squared radius, 1^2+1^2 = 2
  std::cout << "Inner squared radius: " <<
    CGAL::to_double(ma.squared_inner_radius_numerator()) /
    CGAL::to_double(ma.squared_radii_denominator()) << std::endl;

  // get outer squared radius, 2^2+2^2 = 8
  std::cout << "Outer squared radius: " <<
    CGAL::to_double(ma.squared_outer_radius_numerator()) /
    CGAL::to_double(ma.squared_radii_denominator()) << std::endl;

  return 0;

}
