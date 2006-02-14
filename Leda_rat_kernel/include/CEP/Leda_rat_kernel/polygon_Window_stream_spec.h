// provide specialized output operator for rational LEDA kernel traits ...

#ifndef CEP_WINDOW_STREAM_POLYGON_SPEC_2_H
#define CEP_WINDOW_STREAM_POLYGON_SPEC_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/leda_rational.h>
#include <LEDA/rat_window.h>

CGAL_BEGIN_NAMESPACE

template <>
Window_stream&
operator<<(Window_stream& ws,
           const Polygon_2<CGAL::leda_rat_kernel_traits,std::list<leda_rat_point> > &polygon)
{
  typedef Polygon_2<CGAL::leda_rat_kernel_traits,std::list<leda_rat_point> >::Edge_const_circulator EI;
  typedef leda_rat_point  Point_2;

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

CGAL_END_NAMESPACE

#endif 

