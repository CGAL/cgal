// Example program for the barycenter() function for 2D and 3D points.

#include <CGAL/Cartesian.h>
#include <CGAL/barycenter.h>

#include <list>
#include <iostream>
#include <utility>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Point_2           Point_2;
typedef K::Point_3           Point_3;

int main()
{
  // barycenter of 2D weighted points
  std::list<std::pair<Point_2, FT> > points_2;
  points_2.push_back(std::make_pair(Point_2(1.0, 0.0),  1.0));
  points_2.push_back(std::make_pair(Point_2(2.0, 2.0),  2.0));
  points_2.push_back(std::make_pair(Point_2(3.0, 5.0), -2.0));

  Point_2 c2 = CGAL::barycenter(points_2.begin(), points_2.end());
  std::cout << c2 << std::endl;

  // barycenter of 3D weighted points
  std::list<std::pair<Point_3, FT> > points_3;
  points_3.push_back(std::make_pair(Point_3(1.0, 0.0, 0.5),  1.0));
  points_3.push_back(std::make_pair(Point_3(2.0, 2.0, 1.2),  2.0));
  points_3.push_back(std::make_pair(Point_3(3.0, 5.0, 4.5), -5.0));

  Point_3 c3 = CGAL::barycenter(points_3.begin(), points_3.end());
  std::cout << c3 << std::endl;

  return 0;
}
