
#include <CGAL/basic.h>
#include <iostream>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <LEDA/rat_window.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/leda_rational.h>
#include <vector>
#include <list>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_2                                             Point;
typedef std::vector<Point>                                     Container;
typedef CGAL::Polygon_2<K, Container>                          Polygon;

//--------------------------------------------------------------------------//
//                   DrawPolyline
//--------------------------------------------------------------------------//
// draws the polygon 'polygon' on window W, except for the last segment

template <class Traits, class Container>
void DrawPolyline(CGAL::Window_stream& W, const CGAL::Polygon_2<Traits,Container>& polygon)
{
  typedef typename CGAL::Polygon_2<Traits,Container>::Vertex_const_iterator VI;
  typedef typename CGAL::Polygon_2<Traits,Container>::Edge_const_iterator EI;

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

template<class CONTAINER>
CGAL::Window_stream&
draw_fcn(CGAL::Window_stream& ws,
           const CGAL::Polygon_2<CGAL::leda_rat_kernel_traits, CONTAINER>& polygon)
{
  typedef typename CGAL::Polygon_2<CGAL::leda_rat_kernel_traits, CONTAINER>  MY_POLY;
  typedef typename MY_POLY::Edge_const_circulator EI;
  typedef leda_rat_point                 Point_2;
  
  EI e = polygon.edges_circulator();
  
#if defined(CGAL_USE_CGAL_WINDOW)
  CGAL::color cl = ws.get_fill_color();
  
  if (cl != CGAL::invisible) { 
    std::list<CGAL::window_point> LP;
      
    if (e != NULL) {
     EI end = e;
     do {
      Point_2 p = (*e).source();
      double x = CGAL::to_double(p.xcoord());
      double y = CGAL::to_double(p.ycoord());      
     
      LP.push_back(CGAL::window_point(x,y));
      ++e;
      } while (e != end);
    }
    
    ws.draw_filled_polygon(LP,cl);    
  }
#else  
  leda_color cl = ws.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled polygon ...
    leda_list<leda_point> LP;
      
    if (e != NULL) {
     EI end = e;
     do {
      Point_2 p = (*e).source();
      double x = CGAL::to_double(p.xcoord());
      double y = CGAL::to_double(p.ycoord());      
     
      LP.append(leda_point(x,y));
      ++e;
      } while (e != end);
    }
    
    ws.draw_filled_polygon(LP,cl);    
  }
#endif  
  else {
  if (e != NULL) {
    EI end = e;
    do {
      ws << *e;
      ws << (*e).source();
      ++e;
      } while (e != end);
    }
  }
  
  return ws;
}


template <class Traits, class Container>
void DrawPolygon(CGAL::Window_stream& W, const CGAL::Polygon_2<Traits,Container>& polygon)
{
  // operator does not work because of missing x(), y() in the leda_rat_points
  draw_fcn(W,polygon);
}

//--------------------------------------------------------------------------//
//                   PrintPolygonInfo
//--------------------------------------------------------------------------//
// prints some information about the polygon P to cerr

template <class Traits, class Container>
void PrintPolygonInfo(const CGAL::Polygon_2<Traits,Container>& P)
{
  cerr << "Polygon information:" << endl;
  cerr << endl;
  cerr << "  P.size()               = " << P.size() << endl;
  cerr << "  P.is_empty()           = " << int(P.is_empty()) << endl;

  cerr << "  P.is_simple()          = " << int(P.is_simple()) << endl;
  cerr << "  P.is_convex()          = " << int(P.is_convex()) << endl;

  CGAL::Orientation o = P.orientation();
  cerr << "  P.orientation()        = ";
  switch (o) {
    case CGAL::CLOCKWISE       : cerr << "clockwise" << endl; break;
    case CGAL::COUNTERCLOCKWISE: cerr << "counter clockwise" << endl; break;
    case CGAL::COLLINEAR       : cerr << "collinear" << endl; break;
  }

//  cerr << "  P.bbox()               = " << P.bbox() << endl;
  cerr << "  P.area()               = " << P.area() << endl;
  cerr << "  P.left_vertex()        = " << *P.left_vertex() << endl;
  cerr << "  P.right_vertex()       = " << *P.right_vertex() << endl;
  cerr << "  P.top_vertex()         = " << *P.top_vertex() << endl;
  cerr << "  P.bottom_vertex()      = " << *P.bottom_vertex() << endl;
}

//--------------------------------------------------------------------------//
//                   PolygonDemo
//--------------------------------------------------------------------------//

void PolygonDemo(CGAL::Window_stream &W)
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
      W << CGAL::RED;
      DrawPolyline(W, polygon);
    }

    if (right_button_pressed) {
      W.clear();
      W << CGAL::RED;
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
      W << CGAL::GREEN;
      DrawPolygon(W, polygon);
      W << CGAL::RED << p;

      CGAL::Bounded_side bside   = polygon.bounded_side(p);
      switch (bside) {
        case CGAL::ON_BOUNDED_SIDE:
          cout << "  the point is inside the polygon" << endl; break;

        case CGAL::ON_BOUNDARY:
          cout << "  the point is on the boundary of the polygon" << endl; break;

        case CGAL::ON_UNBOUNDED_SIDE:
          cout << "  the point is outside the polygon" << endl; break;
      }
    }

    if (right_button_pressed) {
      W.clear();
      W << CGAL::GREEN;
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
  CGAL::Window_stream W(winx, winy); // physical window size

  double xmin = 0;
  double xmax = 500;
  W.init(xmin, xmax, xmin);         // logical window size

  W.set_show_coordinates(false);
  W.display();

  PolygonDemo(W);

  return 0;
}

