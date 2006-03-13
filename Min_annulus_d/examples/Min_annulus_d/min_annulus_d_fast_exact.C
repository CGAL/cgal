// file: examples/Min_annulus_d/min_annulus_d_fast_exact.C

// computes the smallest enclosing annulus of two point
// sets on nested squares in R^2, using double as input
// type and CGAL::Double as exact internal type; this
// is guaranteed to have no roundoff errors. Note: CGAL::Double
// is based on GMP but not yet an official CGAL number type; in 
// this respect, the example represents experimental code 
#include <CGAL/QP_solver/gmp_double.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/QP_solver/Double.h>
#include <CGAL/Min_annulus_d.h>
#include <CGAL/Optimisation_d_traits_2.h>

typedef CGAL::Simple_cartesian<double>    K;  
typedef K::Point_2                        Point;
typedef CGAL::Optimisation_d_traits_2<K, CGAL::Double, double>  
                                          Traits;
typedef CGAL::Min_annulus_d<Traits>       Min_annulus;


int main()
{
  // points on the squares [-1,1]^2 and [-2,2]^2
  Point P[8] = { Point(-1,-1), Point(-1,1), Point(1,-1), Point(1,1),
                 Point(-2,-2), Point(-2,2), Point(2,-2), Point(2,2)};

  Min_annulus ma(P, P+8); 

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
