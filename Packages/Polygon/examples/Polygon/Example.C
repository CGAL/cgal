//-----------------------------------------------------------------------//
// This is just a simple example that demonstrates how to use the
// class CGAL_Polygon_2.
//-----------------------------------------------------------------------//

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>

#include <list.h>

typedef CGAL_Cartesian<double> R;
typedef CGAL_Polygon_traits_2<R> Traits;
typedef CGAL_Point_2<R> Point;
typedef CGAL_Polygon_2<Traits, list<Point> > Polygon;
typedef CGAL_Polygon_2<Traits, list<Point> >::Vertex_iterator VertexIterator;
typedef CGAL_Polygon_2<Traits, list<Point> >::Edge_const_iterator EdgeIterator;

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  // create a polygon and put some points in it
  Polygon p;
  p.push_back(Point(0,0));
  p.push_back(Point(4,0));
  p.push_back(Point(4,4));
  p.push_back(Point(2,2));
  p.push_back(Point(0,4));

  CGAL_set_pretty_mode(cout);
  cout << "created the polygon p:" << endl;
  cout << p << endl;
  cout << endl;

  // determine some properties of the polygon
  bool IsSimple    = p.is_simple();
  bool IsConvex    = p.is_convex();
  bool IsClockwise = (p.orientation() == CGAL_CLOCKWISE);
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
  Point q(1,1);
  cout << "created point q = " << q << endl;
  cout << endl;

  bool IsInside = (p.bounded_side(q) == CGAL_ON_BOUNDED_SIDE);
  cout << "point q is";
  if (!IsInside) cout << " not";
  cout << " inside polygon p." << endl;
  cout << endl;

  // traverse the vertices and the edges
  int n=0;
  for (VertexIterator vi = p.vertices_begin(); vi != p.vertices_end(); ++vi)
    cout << "vertex " << n++ << " = " << *vi << endl;
  cout << endl;

  n=0;
  for (EdgeIterator ei = p.edges_begin(); ei != p.edges_end(); ++ei)
    cout << "edge " << n++ << " = " << *ei << endl;

  return 0;
}

