#include <CGAL/Simple_cartesian.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/iterator.h>
#include <iostream>
#include <algorithm>
#include <array>


typedef double                                      FT;

typedef CGAL::Simple_cartesian<FT>                  Kernel;

typedef Kernel::Point_2                             Point;
typedef std::array<Point,4>                          Cont;

int main()
{
  int n = 10;
  int p = 3;

  Cont points = { Point(0, -1), Point(1, -2),  Point(3,2),  Point(4,1) };

  FT p_radius;

  CGAL::rectangular_p_center_2(
    points.begin(), points.end(), CGAL::Emptyset_iterator(), p_radius, p);
  std::cout << "\n\n" << p << "-radius = " << p_radius << std::endl;

  return 0;
}
