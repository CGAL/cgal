// file: examples/Min_annulus_d/min_annulus_d.C

// computes the smallest enclosing annulus of two point
// sets on nested squares in R^2, using double as input
// type and internal type (prone to roundoff errors)

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_annulus_d.h>
#include <CGAL/Optimisation_d_traits_2.h>

typedef CGAL::Simple_cartesian<double>    K;  
typedef K::Point_2                        Point;
typedef CGAL::Optimisation_d_traits_2<K>  Traits;
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
    std::cout << " " << *coord_it;
  std::cout << std::endl;

  // get inner squared radius, 1^2+1^2 = 2
  std::cout << "Inner squared radius: " <<
    ma.squared_inner_radius_numerator() /
    ma.squared_radii_denominator() << std::endl;

  // get outer squared radius, 2^2+2^2 = 8
  std::cout << "Outer squared radius: " <<
    ma.squared_outer_radius_numerator() /
    ma.squared_radii_denominator() << std::endl;

  return 0;

} 
