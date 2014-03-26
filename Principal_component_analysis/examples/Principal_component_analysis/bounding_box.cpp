// Example program for the bounding_box() function for 2D points and 3D points.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>

#include <vector>
#include <iostream>

typedef double                     FT;
typedef CGAL::Simple_cartesian<FT> K;
typedef K::Point_2                 Point_2;
typedef K::Point_3                 Point_3;

int main()
{
  // axis-aligned bounding box of 2D points
  std::vector<Point_2> points_2;
  points_2.push_back(Point_2(1.0, 0.0));
  points_2.push_back(Point_2(2.0, 2.0));
  points_2.push_back(Point_2(3.0, 5.0));

  K::Iso_rectangle_2 c2 = CGAL::bounding_box(points_2.begin(), points_2.end());
  std::cout << c2 << std::endl;

  // axis-aligned bounding box of 3D points
  std::vector<Point_3> points_3;
  points_3.push_back(Point_3(1.0, 0.0, 0.5));
  points_3.push_back(Point_3(2.0, 2.0, 1.2));
  points_3.push_back(Point_3(3.0, 5.0, 4.5));

  K::Iso_cuboid_3 c3 = CGAL::bounding_box(points_3.begin(), points_3.end());
  std::cout << c3 << std::endl;

  return 0;
}
