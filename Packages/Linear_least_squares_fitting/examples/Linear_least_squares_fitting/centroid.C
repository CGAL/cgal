// Example program for the centroid() function for 2D and 3D points.

#include <CGAL/Cartesian.h>
#include <CGAL/centroid.h>

#include <list>
#include <iostream>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Point_2           Point_2;
typedef K::Point_3           Point_3;

int main()
{
  // centroid of 2D points
  std::list<Point_2> points_2;
  points_2.push_back(Point_2(1.0, 0.0));
  points_2.push_back(Point_2(2.0, 2.0));
  points_2.push_back(Point_2(3.0, 5.0));

  Point_2 c2 = CGAL::centroid(points_2.begin(), points_2.end());
  std::cout << c2 << std::endl;

  // centroid of 3D points
  std::list<Point_3> points_3;
  points_3.push_back(Point_3(1.0, 0.0, 0.5));
  points_3.push_back(Point_3(2.0, 2.0, 1.2));
  points_3.push_back(Point_3(3.0, 5.0, 4.5));

  Point_3 c3 = CGAL::centroid(points_3.begin(), points_3.end());
  std::cout << c3 << std::endl;

  return 0;
}
