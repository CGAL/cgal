// computes the distance between two cubes in R^3 using double
// as input type and some internal EXACT floating point type
#include <iostream>
#include <cassert>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_3.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

// use an EXACT kernel...
typedef CGAL::Homogeneous<ET>                 K;
typedef K::Point_3                            Point;
// ...and the traits class based on the exact kernel
typedef CGAL::Polytope_distance_d_traits_3<K> Traits;
typedef CGAL::Polytope_distance_d<Traits>     Polytope_distance;


int main()
{
  // the cube [0,1]^3
  Point P[8] = { Point(0,0,0), Point(0,0,1), Point(0,1,0), Point(0,1,1),
                 Point(1,0,0), Point(1,0,1), Point(1,1,0), Point(1,1,1)};

  // the cube [2,3]^3
  Point Q[8] = { Point(2,2,2), Point(2,2,3), Point(2,3,2), Point(2,3,3),
                 Point(3,2,2), Point(3,2,3), Point(3,3,2), Point(3,3,3)};

  Polytope_distance pd(P, P+8, Q, Q+8);
  assert (pd.is_valid());

  // get squared distance (2,2,2)-(1,1,1))^2 = 3
  std::cout << "Squared distance: " <<
    CGAL::to_double (pd.squared_distance_numerator()) /
    CGAL::to_double (pd.squared_distance_denominator()) << std::endl;

  // get points that realize the distance
  Polytope_distance::Coordinate_iterator  coord_it;

  std::cout << "p:"; // homogeneous point from first cube, (1,1,1,1)
  for (coord_it = pd.realizing_point_p_coordinates_begin();
       coord_it != pd.realizing_point_p_coordinates_end();
       ++coord_it)
    std::cout << " " << *coord_it;
  std::cout << std::endl;

  std::cout << "q:"; // homogeneous point from second cube, (2,2,2,1)
  for (coord_it = pd.realizing_point_q_coordinates_begin();
       coord_it != pd.realizing_point_q_coordinates_end();
       ++coord_it)
    std::cout << " " << *coord_it;
  std::cout << std::endl;

  return 0;

}
