//-----------------------------------------------------------------------//
// This is a small graphical demo program for polygons.
// N.B. The program sometimes requires keyboard input!
//-----------------------------------------------------------------------//

#include <CGAL/basic.h>
#include <fstream.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <vector.h>

typedef CGAL_Cartesian<double> R;
typedef CGAL_Point_2<R> Point;
typedef CGAL_Polygon_traits_2<R> Traits;
typedef vector<Point> Container;
typedef CGAL_Polygon_2<Traits, Container> Polygon;

//--------------------------------------------------------------------------//
//                   DrawPolyline
//--------------------------------------------------------------------------//
// draws the polygon 'polygon' on window W, except for the last segment

template <class Traits, class Container>
void DrawPolyline(CGAL_Window_stream& W, const CGAL_Polygon_2<Traits,Container>& polygon)
{
  typedef typename CGAL_Polygon_2<Traits,Container>::Vertex_const_iterator VI;
  typedef typename CGAL_Polygon_2<Traits,Container>::Edge_const_iterator EI;

  VI v_begin = polygon.vertices_begin();
  VI v_end   = polygon.vertices_end();
  for (VI v = v_begin; v != v_end; ++v)
    W << *v;

  EI e_begin = polygon.edges_begin();
  EI e_end   = polygon.edges_end();
  if (e_end != e_begin) --e_end; // don't draw the last edge
  for (EI e = e_begin; e != e_end; ++e)
    W << *e;
}

//--------------------------------------------------------------------------//
//                   DrawPolygon
//--------------------------------------------------------------------------//
// draws the polygon 'polygon' on window W

template <class Traits, class Container>
void DrawPolygon(CGAL_Window_stream& W, const CGAL_Polygon_2<Traits,Container>& polygon)
{
  W << polygon;
}

//--------------------------------------------------------------------------//
//                   PrintPolygonInfo
//--------------------------------------------------------------------------//
// prints some information about the polygon P to cerr

template <class Traits, class Container>
void PrintPolygonInfo(const CGAL_Polygon_2<Traits,Container>& P)
{
  cerr << "Polygon information:" << endl;
  cerr << endl;
  cerr << "  P.size()               = " << P.size() << endl;
  cerr << "  P.is_empty()           = " << int(P.is_empty()) << endl;

  cerr << "  P.is_simple()          = " << int(P.is_simple()) << endl;
  cerr << "  P.is_convex()          = " << int(P.is_convex()) << endl;

  CGAL_Orientation o = P.orientation();
  cerr << "  P.orientation()        = ";
  switch (o) {
    case CGAL_CLOCKWISE       : cerr << "clockwise" << endl; break;
    case CGAL_COUNTERCLOCKWISE: cerr << "counter clockwise" << endl; break;
    case CGAL_COLLINEAR       : cerr << "collinear" << endl; break;
  }

  cerr << "  P.bbox()               = " << P.bbox() << endl;
  cerr << "  P.area()               = " << P.area() << endl;
  cerr << "  P.left_vertex()        = " << *P.left_vertex() << endl;
  cerr << "  P.right_vertex()       = " << *P.right_vertex() << endl;
  cerr << "  P.top_vertex()         = " << *P.top_vertex() << endl;
  cerr << "  P.bottom_vertex()      = " << *P.bottom_vertex() << endl;
}

//--------------------------------------------------------------------------//
//                   PolygonDemo
//--------------------------------------------------------------------------//

void PolygonDemo(CGAL_Window_stream &W)
{
  cerr << "Enter points with the left button" << endl;
  cerr << "Right button terminates input of points" << endl;

  Polygon polygon;

  while(1) {
    double x, y;
    int b = W.get_mouse(x,y);
    bool left_button_pressed  = (b == MOUSE_BUTTON(1));
    bool right_button_pressed = (b == MOUSE_BUTTON(3));

    Point p(x,y);

    if (left_button_pressed) {
      polygon.push_back(p);
      W.clear();
      W << CGAL_RED;
      DrawPolyline(W, polygon);
    }

    if (right_button_pressed) {
      W.clear();
      W << CGAL_RED;
      DrawPolygon(W, polygon);
      break;
    }
  }

  PrintPolygonInfo(polygon);

  cerr << "Enter a point with the left button" << endl;
  cerr << "Right button terminates input of points" << endl;
  while (1) {
    double x, y;
    int b = W.get_mouse(x,y);
    bool left_button_pressed  = (b == MOUSE_BUTTON(1));
    bool right_button_pressed = (b == MOUSE_BUTTON(3));

    Point p(x,y);

    if (left_button_pressed) {
      W.clear();
      W << CGAL_GREEN;
      DrawPolygon(W, polygon);
      W << CGAL_RED << p;

      CGAL_Bounded_side bside   = polygon.bounded_side(p);
      switch (bside) {
        case CGAL_ON_BOUNDED_SIDE:
          cout << "  the point is inside the polygon" << endl; break;

        case CGAL_ON_BOUNDARY:
          cout << "  the point is on the boundary of the polygon" << endl; break;

        case CGAL_ON_UNBOUNDED_SIDE:
          cout << "  the point is outside the polygon" << endl; break;
      }
    }

    if (right_button_pressed) {
      W.clear();
      W << CGAL_GREEN;
      DrawPolygon(W, polygon);
      break;
    }
  }

  cerr << "press enter to quit ..." << endl;
  cin.get();
}

int main()
{
  int winx = 500;
  int winy = 500;
  CGAL_Window_stream W(winx, winy); // physical window size

  double xmin = 0;
  double xmax = 500;
  W.init(xmin, xmax, xmin);         // logical window size

  W.set_show_coordinates(false);

  PolygonDemo(W);

  return 0;
}

