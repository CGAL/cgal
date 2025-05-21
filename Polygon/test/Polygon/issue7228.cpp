#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <vector>
#include <array>
#include <cassert>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef Polygon_2::Vertex_circulator Vertex_circulator;

int main()
{
  std::array<Point,4> points =  { Point(0,0), Point(1,0), Point(1,1), Point(0,1) };
  Polygon_2 poly(points.begin(), points.end());

  Vertex_circulator vc = poly.vertices_circulator();

  ++vc;
  ++vc;
  ++vc;

  vc = poly.erase(vc);

  assert(*vc == Point(0,0));

  return 0;
}
