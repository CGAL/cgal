// Example program for the bounding_box() function for 2D points,
// 3D points and 3D segments.

#include <CGAL/Cartesian.h>
#include <CGAL/bounding_box.h>

#include <list>
#include <iostream>

typedef double              FT;
typedef CGAL::Cartesian<FT> K;
typedef K::Point_2          Point_2;
typedef K::Point_3          Point_3;
typedef K::Point_3          Point_3;
typedef K::Segment_3        Segment_3;

int main()
{
  // axis-aligned bounding box of 2D points
  std::list<Point_2> points_2;
  points_2.push_back(Point_2(1.0, 0.0));
  points_2.push_back(Point_2(2.0, 2.0));
  points_2.push_back(Point_2(3.0, 5.0));

  K::Iso_rectangle_2 c2 = CGAL::bounding_box(points_2.begin(), points_2.end());
  std::cout << c2 << std::endl;

  // axis-aligned bounding box of 3D points
  std::list<Point_3> points_3;
  points_3.push_back(Point_3(1.0, 0.0, 0.5));
  points_3.push_back(Point_3(2.0, 2.0, 1.2));
  points_3.push_back(Point_3(3.0, 5.0, 4.5));

  K::Iso_cuboid_3 c3 = CGAL::bounding_box(points_3.begin(), points_3.end());
  std::cout << c3 << std::endl;

  // axis-aligned bounding box of 3D segments
  std::list<Segment_3> segments_3;
  Point_3 p(1.0, 2.0, 3.0);
  Point_3 q(4.0, 5.0, 6.0);
  Point_3 r(3.0, 3.0, 0.5);
  segments_3.push_back(Segment_3(p,q));
  segments_3.push_back(Segment_3(p,r));

  c3 = CGAL::bounding_box(segments_3.begin(), segments_3.end());
  std::cout << c3 << std::endl;

  return 0;
}
