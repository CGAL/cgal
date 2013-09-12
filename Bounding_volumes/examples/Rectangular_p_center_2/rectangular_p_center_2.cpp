#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/algorithm.h>
#include <iostream>
#include <algorithm>
#include <vector>

typedef double                                      FT;

typedef CGAL::Simple_cartesian<FT>                  Kernel;

typedef Kernel::Point_2                             Point;
typedef std::vector<Point>                          Cont;
typedef CGAL::Random_points_in_square_2<Point>      Generator;
typedef CGAL::Ostream_iterator<Point,std::ostream>  OIterator;

int main()
{
  int n = 10;
  int p = 2;
  OIterator cout_ip(std::cout);
  CGAL::set_pretty_mode(std::cout);

  Cont points;
  CGAL::cpp11::copy_n(Generator(1), n, std::back_inserter(points));
  std::cout << "Generated Point Set:\n";
  std::copy(points.begin(), points.end(), cout_ip);

  FT p_radius;
  std::cout << "\n\n" << p << "-centers:\n";
  CGAL::rectangular_p_center_2(
    points.begin(), points.end(), cout_ip, p_radius, 3);
  std::cout << "\n\n" << p << "-radius = " << p_radius << std::endl;

  return 0;
}
