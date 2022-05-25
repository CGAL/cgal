#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

#include <vector>
#include <iostream>
#include <cassert>
#include <iterator>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;


typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::General_polygon_with_holes_2<Polygon_2> Polygon_with_holes_2;

int main()
{
  std::array<Point,4> outer = { Point(0, 0), Point(10, 0), Point(10, 10), Point(0, 10) };
  std::array<Point, 4> hole1 = { Point(1, 1), Point(1, 2), Point(2, 2), Point(2, 1) };
  std::array<Point, 4> hole2 = { Point(3, 3), Point(3, 4), Point(4, 4), Point(4, 3) };

  std::vector<Polygon_2> holes;
  holes.reserve(2);
  holes.emplace_back(hole1.begin(), hole1.end());
  holes.emplace_back(hole2.begin(), hole2.end());

  Polygon_2 pouter(outer.begin(), outer.end());

  Polygon_with_holes_2 pwh(std::move(pouter), std::move_iterator<std::vector<Polygon_2>::iterator>(holes.begin()), std::move_iterator<std::vector<Polygon_2>::iterator>(holes.end()));

  std::cout << pwh << std::endl;
  assert(pouter.is_empty());
  assert(holes[0].is_empty());
  assert(holes[1].is_empty());
  return 0;
}
