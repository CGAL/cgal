//-----------------------------------------------------------------------//
// This example demonstrates how to define a polygon with a custom point
// type with additional attributes.
//
// For this the following is needed:
//
// 1) define a point, vector and segment type
// 2) define a polygon traits class using these types
//
// In this example the type MyPoint derives from CGAL::Point_2. As a result
// most of the work from 1) and 2) isn't necessary. The types CGAL::Segment_2
// and CGAL::Vector_2 and an existing traits class can be reused.
//-----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>

enum MyColor { Red, Green, Blue };

using namespace CGAL;

template <class R>
class MyPoint: public Point_2<R>
{
  public:
    MyColor color; // additional attribute for points

    MyPoint(): color(Red) {}
    MyPoint(const typename R::RT &x, const typename R::RT &y, MyColor c)
      : Point_2<R>(x,y), color(c) {}
};
// N.B. The class MyPoint derives from CGAL::Point_2. This is not always a
// good idea, since the class CGAL::Point_2 doesn't have a virtual destructor.

typedef Cartesian<double> R;
typedef MyPoint<R> Point;

#include <CGAL/Polygon_2.h>
#include <list>

typedef Polygon_traits_2_aux<R, R::FT, Point> Traits;
// The class MyPoint derives from CGAL::Point, so the polygon traits class
// CGAL::Polygon_traits_2_aux can be reused.

typedef Polygon_2<Traits, std::list<Point> > Polygon;
typedef Polygon_2<Traits, std::list<Point> >::Vertex_iterator VI;
typedef Polygon_2<Traits, std::list<Point> >::Edge_const_iterator EI;

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  // create a polygon and put some points in it
  Polygon p;
  p.push_back(Point(0,0,Red));
  p.push_back(Point(4,0,Blue));
  p.push_back(Point(4,4,Green));
  p.push_back(Point(2,2,Red));
  p.push_back(Point(0,4,Red));

  CGAL::set_pretty_mode(std::cout);
  std::cout << "created the polygon p:" << std::endl;
  std::cout << p << std::endl;
  std::cout << std::endl;

  // determine some properties of the polygon
  bool IsSimple    = p.is_simple();
  bool IsConvex    = p.is_convex();
  bool IsClockwise = (p.orientation() == CGAL::CLOCKWISE);
  double Area      = p.area();

  std::cout << "polygon p is";
  if (!IsSimple) std::cout << " not";
  std::cout << " simple." << std::endl;

  std::cout << "polygon p is";
  if (!IsConvex) std::cout << " not";
  std::cout << " convex." << std::endl;

  std::cout << "polygon p is";
  if (!IsClockwise) std::cout << " not";
  std::cout << " clockwise oriented." << std::endl;

  std::cout << "the area of polygon p is " << Area << std::endl;
  std::cout << std::endl;

  // apply some algorithms
  Point q(1,1,Blue);
  std::cout << "created point q = " << q << std::endl;
  std::cout << std::endl;

  bool IsInside = (p.bounded_side(q) == CGAL::ON_BOUNDED_SIDE);
  std::cout << "point q is";
  if (!IsInside) std::cout << " not";
  std::cout << " inside polygon p." << std::endl;
  std::cout << std::endl;

  // traverse the vertices and the edges
  int n=0;
  for (VI vi = p.vertices_begin(); vi != p.vertices_end(); ++vi)
    std::cout << "vertex " << n++ << " = " << *vi << std::endl;
  std::cout << std::endl;

  n=0;
  for (EI ei = p.edges_begin(); ei != p.edges_end(); ++ei)
    std::cout << "edge " << n++ << " = " << *ei << std::endl;

  return 0;
}

