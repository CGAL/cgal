#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <iostream>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>           Traits;
typedef Traits::FT                                                              NT;
typedef Traits::Point_2                                                         Point;
typedef Traits::Circle_2                                                        Circle;
typedef Traits::Euclidean_line_2                                                Line;
typedef Traits::Construct_intersection_2                                        Construct_intersection_2;

int main(int /*argc*/, char** /*argv*/)
{
  NT F2(2);

  Line ell1(NT(1),  sqrt(F2 - sqrt(F2)), NT(0));
  Line ell2(NT(1), -sqrt(F2 + sqrt(F2)), NT(0));

  // In the Poincaré disk, all Euclidean lines pass through the origin
  NT sx = NT(0);
  NT sy = NT(0);
  Point sol(sx, sy);

  Point p = Construct_intersection_2()(ell1, ell2);
  std::cout << "line-line intersection: " << p << ", solution is: " << sol << std::endl;
  assert(p == sol);
  std::cout << "The solution is exact!" << std::endl;

  Point pc1(NT(2)/NT(10), NT(5)/NT(10));
  Point pc2(-NT(3)/NT(10), NT(15)/NT(100));
  Point pc3(NT(20)/NT(29), NT(50)/NT(29));
  Circle c1(pc1, pc2, pc3);
  std::pair<Point, Point> ipt = Construct_intersection_2()(ell1, c1);
  Point ip1 = ipt.first, ip2 = ipt.second;

  NT sq(sqrt(NT(2)));
  NT nsq( sqrt( NT(2) - sqrt(NT(2)) ) );
  NT sol1x( nsq*(NT(2438)+NT(1451)*nsq-sqrt(NT(3933846)-NT(31801)*sq+NT(7075076)*nsq))/(-NT(4320)+NT(1440)*sq) );
  NT sol1y( (-NT(2438)-NT(1451)*nsq+sqrt(NT(3933846)-NT(31801)*sq+NT(7075076)*nsq))/(-NT(4320)+NT(1440)*sq) );
  Point solution1(sol1x, sol1y);

  std::cout << "Intersection of circle and line: " << ip1 << " and " << ip2 << std::endl;
  std::cout << "Solution: " << solution1 << std::endl;
  assert(ip1 == solution1 || ip2 == solution1);
  std::cout << "The solution is exact!" << std::endl;

  Point pc4(NT(2)/NT(10), -NT(15)/NT(100));
  Circle c2(pc1, pc4, pc3);

  std::pair<Point, Point> cpt = Construct_intersection_2()(c1, c2);
  ip1 = cpt.first;
  ip2 = cpt.second;
  std::cout << "Intersection of circles: " << ip1 << " and " << ip2 << std::endl;
  assert(ip1 == pc1 || ip2 == pc1);
  std::cout << "The solution is exact!" << std::endl;

  return EXIT_SUCCESS;
}
