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

//#ifdef CGAL_CFG_NO_TEMPLATE_FUNCTION_MATCHING
// This is a workaround for the g++ 2.7.2 compiler.
//#include "template_function_matching_workaround.h"
//#endif // CGAL_CFG_NO_TEMPLATE_FUNCTION_MATCHING

#include <CGAL/Polygon_2.h>
#include <list.h>

typedef Polygon_traits_2_aux<R, R::FT, Point> Traits;
// The class MyPoint derives from CGAL::Point, so the polygon traits class
// CGAL::Polygon_traits_2_aux can be reused.

typedef Polygon_2<Traits, list<Point> > Polygon;
typedef Polygon_2<Traits, list<Point> >::Vertex_iterator VI;
typedef Polygon_2<Traits, list<Point> >::Edge_const_iterator EI;

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

  set_pretty_mode(cout);
  cout << "created the polygon p:" << endl;
  cout << p << endl;
  cout << endl;

  // determine some properties of the polygon
  bool IsSimple    = p.is_simple();
  bool IsConvex    = p.is_convex();
  bool IsClockwise = (p.orientation() == CGAL::CLOCKWISE);
  double Area      = p.area();

  cout << "polygon p is";
  if (!IsSimple) cout << " not";
  cout << " simple." << endl;

  cout << "polygon p is";
  if (!IsConvex) cout << " not";
  cout << " convex." << endl;

  cout << "polygon p is";
  if (!IsClockwise) cout << " not";
  cout << " clockwise oriented." << endl;

  cout << "the area of polygon p is " << Area << endl;
  cout << endl;

  // apply some algorithms
  Point q(1,1,Blue);
  cout << "created point q = " << q << endl;
  cout << endl;

  bool IsInside = (p.bounded_side(q) == CGAL::ON_BOUNDED_SIDE);
  cout << "point q is";
  if (!IsInside) cout << " not";
  cout << " inside polygon p." << endl;
  cout << endl;

  // traverse the vertices and the edges
  int n=0;
  for (VI vi = p.vertices_begin(); vi != p.vertices_end(); ++vi)
    cout << "vertex " << n++ << " = " << *vi << endl;
  cout << endl;

  n=0;
  for (EI ei = p.edges_begin(); ei != p.edges_end(); ++ei)
    cout << "edge " << n++ << " = " << *ei << endl;

  return 0;
}

