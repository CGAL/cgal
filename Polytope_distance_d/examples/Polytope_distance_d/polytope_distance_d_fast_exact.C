// file: examples/Polytope_distance_d/polytope_distance_fast_exact.C

// computes the distance between two cubes in R^3 using double
// as input type and CGAL::Double as exact internal type; this
// is guaranteed to have no roundoff errors. Note: CGAL::Double
// is based on GMP but not yet an official CGAL number type; in 
// this respect, the example represents experimental code
#include <iostream>
#include <CGAL/QP_solver/gmp_double.h> // will become CGAL number type
#include <CGAL/Simple_cartesian.h>
#include <CGAL/QP_solver/Double.h>     // will become CGAL number type
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Optimisation_d_traits_3.h>

typedef CGAL::Simple_cartesian<double>    K;  
typedef K::Point_3                        Point;
typedef CGAL::Optimisation_d_traits_3<K, CGAL::Double, double>  
                                          Traits;
typedef CGAL::Polytope_distance_d<Traits> Polytope_distance;

int main()
{
  // the cube [0,1]^3
  Point P[8] = { Point(0,0,0), Point(0,0,1), Point(0,1,0), Point(0,1,1),
                 Point(1,0,0), Point(1,0,1), Point(1,1,0), Point(1,1,1)};

  // the cube [2,3]^3
  Point Q[8] = { Point(2,2,2), Point(2,2,3), Point(2,3,2), Point(2,3,3),
                 Point(3,2,2), Point(3,2,3), Point(3,3,2), Point(3,3,3)};

  Polytope_distance pd(P, P+8, Q, Q+8); 

  // get squared distance (2,2,2)-(1,1,1))^2 = 3
  std::cout << "Squared Distance: " <<
    CGAL::to_double(pd.squared_distance_numerator()) /
    CGAL::to_double(pd.squared_distance_denominator()) << std::endl;

  // get points that realize the distance
  Polytope_distance::Coordinate_iterator  coord_it;

  std::cout << "p:"; // homogeneous point from first cube, (1,1,1,1)
  for (coord_it = pd.realizing_point_p_coordinates_begin();
       coord_it != pd.realizing_point_p_coordinates_end();
       ++coord_it)
    std::cout << " " <<  CGAL::to_double(*coord_it);
  std::cout << std::endl;

  std::cout << "q:"; // homogeneous point from second cube, (2,2,2,1)
  for (coord_it = pd.realizing_point_q_coordinates_begin();
       coord_it != pd.realizing_point_q_coordinates_end();
       ++coord_it)
    std::cout << " " <<  CGAL::to_double(*coord_it);
  std::cout << std::endl;

  return 0;
}
  
