#include <CGAL/Simple_cartesian.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Line_2 Line_2;

int main()
{
  Line_2 line(-4.2885603045067812e-18, 1, 250.73609999999996);

  Point_2 point(35.306000000000004, 250.69800000000001);
  std::cout.precision(17);
  std::cout << line.projection(point) << std::endl;

  return 0;
}
