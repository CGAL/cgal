//  -*- Mode: c++ -*-
// ============================================================================
// 
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-1.0 $
// release_date  : $CGAL_Date: 1998/09/12 $
//
// file          : demo/BooleanOperations/include/CGAL/IO-old/Window_stream.h
// source        : demo/BooleanOperations/include/CGAL/IO-old/Window_stream.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                        Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#ifndef CGAL_WINDOW_STREAM_H
#define CGAL_WINDOW_STREAM_H


#include <CGAL/IO/Color.h>

#include <LEDA/window.h>
#include <LEDA/REDEFINE_NAMES.h>
typedef window LEDA_Window;
#include <LEDA/UNDEFINE_NAMES.h>

class Window_stream : public LEDA_Window  {

public:
  Window_stream(int xpix = 300,
                     int ypix = 300,
                     int xpos = 400,
                     int ypos = 400)
    : LEDA_Window(xpix, ypix)
  {
    set_frame_label("Window_stream");
    LEDA_Window::display(xpos, ypos);
  }



  void init(const Bbox_2 &b = Bbox_2(0, 0, 100, 100))
  {
    LEDA_Window::init(b.xmin(), b.xmax(), b.ymin());
    set_fg_color(black);
  }

  void init(double xmin, double xmax, double ymin)
  {
    LEDA_Window::init(xmin, xmax, ymin);
    set_fg_color(black);
  }

  int read_mouse(double& x, double& y)
  {
    button = LEDA_Window::read_mouse(x,y);
    return button;
  }

  int last_button_pressed()
  {
    return button;
  }
  // we keep the foreground color in a field, as LEDA's window
  // can't store it in the shared library version
private:
  int button;

};


#ifdef CGAL_POINT_2_H
template< class R >
Window_stream& operator<<(Window_stream &w, const Point_2<R> &p)
{
  w.draw_point(to_double(p.x()), to_double(p.y()));
  return w;
}


template< class R >
Window_stream& operator>>(Window_stream &w, Point_2<R> &p)
{
  double x, y;
  w.read_mouse(x,y);
  w.draw_point(x,y);
  p = Point_2<R>(x,y);
  return w;
}
#endif //  CGAL_POINT_2_H

#ifdef CGAL_LINE_2_H
template< class R >
Window_stream& operator<<(Window_stream &w,const Segment_2<R> &s)
{
  w.draw_segment(to_double(s.source().x()),
                 to_double(s.source().y()),
                 to_double(s.target().x()),
                 to_double(s.target().y()));
  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w,
                               Segment_2<R> &s)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  w.draw_segment(x1,y1, x2, y2);
  s = Segment_2<R>(Point_2<R>(x1,y1),
                        Point_2<R>(x2,y2));
  return w;
}


template< class R >
Window_stream& operator<<(Window_stream &w, const Line_2<R> &l)
{
  Point_2<R> p1 = l.point(),
                  p2 = p1 + l.direction().vector();

  w.draw_line(to_double(p1.x()), to_double(p1.y()),
              to_double(p2.x()), to_double(p2.y()));
  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w, Line_2<R> &l)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  w.draw_line(x1,y1, x2, y2);

  l = Line_2<R>(Point_2<R>(x1,y1),
                     Point_2<R>(x2,y2));
  return w;
}
#endif //CGAL_LINE_2_H

inline Window_stream& operator<<(Window_stream &w,
                                      const Color& c)
{
  w.set_fg_color(color(c.red(), c.green(), c.blue()));

  return w;
}


#ifdef CGAL_ISO_RECTANGLE_2_H
template< class R >
Window_stream& operator<<(Window_stream &w,
                               const Iso_rectangle_2<R> &r)
{
  double xmin = to_double(r.min().x()),
         ymin = to_double(r.min().y()),
         xmax = to_double(r.max().x()),
         ymax = to_double(r.max().y());


  w.draw_segment(xmin, ymin, xmax, ymin);
  w.draw_segment(xmax, ymin, xmax, ymax);
  w.draw_segment(xmax, ymax, xmin, ymax);
  w.draw_segment(xmin, ymax, xmin, ymin);

  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w,
                               Iso_rectangle_2<R> &r)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_rect(x1,y1, x2, y2);
  r = Iso_rectangle_2<R>(Point_2<R>(x1,y1),
                              Point_2<R>(x2,y2));
  w << r;
  return w;
}

#endif // CGAL_ISO_RECTANGLE_2_H

#ifdef CGAL_TRIANGLE_2_H
template< class R >
Window_stream& operator<<(Window_stream &w,
                               const Triangle_2<R> &t)
{
  double x0 = to_double(t.vertex(0).x()),
         y0 = to_double(t.vertex(0).y()),
         x1 = to_double(t.vertex(1).x()),
         y1 = to_double(t.vertex(1).y()),
         x2 = to_double(t.vertex(2).x()),
         y2 = to_double(t.vertex(2).y());


  w.draw_segment(x0, y0, x1, y1);
  w.draw_segment(x1, y1, x2, y2);
  w.draw_segment(x2, y2, x0, y0);

  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w,
                               Triangle_2<R> &t)
{
  double x0, y0, x1, y1, x2, y2;
  w.read_mouse(x0,y0);
  w.read_mouse_seg(x0,y0, x1, y1);
  w.draw_seg(x0,y0, x1, y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  w.draw_seg(x1,y1, x2, y2);
  w.draw_seg(x2,y2, x0, y0);
  r = Triangle_2<R>(Point_2<R>(x0,y0),
                         Point_2<R>(x1,y1),
                         Point_2<R>(x2,y2));
  return w;
}

#endif // CGAL_TRIANGLE_2_H



#if 0

#ifdef CGAL_RAY_2_H

#include <CGAL/Object.h>
#include <CGAL/intersection_2_1.h>

template< class R >
Window_stream& operator<<(Window_stream &w, const Ray_2<R> &r)
{
  Point_2<R> p;
  Segment_2<R> s;
  Point_2<R> point[4];
  Point_2<R> ipoint[2];

  point[0] = Point_2<R>(w.xmin(), w.ymin());
  point[1] = Point_2<R>(w.xmin(), w.ymax());
  point[2] = Point_2<R>(w.xmax(), w.ymax());
  point[3] = Point_2<R>(w.xmax(), w.ymin());

  int no_of_intersections = 0;
  for(int i=0; i < 4; i++) {
    Segment_2<R> seg(point[i],point[(i+1)%4]);
    Object o = intersection(seg, r);
    if ( assign(p, o) ) {
      ipoint[no_of_intersections] = p;
      no_of_intersections++;
    }else if( assign(s, o) ) {
      w.draw_segment(to_double(s.source().x()),
                     to_double(s.source().y()),
                     to_double(s.target().x()),
                     to_double(s.target().y()));
      return w;
    }
  }

  if( no_of_intersections == 1 ) {
    // the start point of the ray is inside the window

    w.draw_segment(to_double(r.source().x()),
                   to_double(r.source().y()),
                   to_double(ipoint[0].x()),
                   to_double(ipoint[0].y()));
  } else if( no_of_intersections == 2 ) {
    w.draw_segment(to_double(ipoint[0].x()),
                   to_double(ipoint[0].y()),
                   to_double(ipoint[1].x()),
                   to_double(ipoint[1].y()));
  }

  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w, Ray_2<R> &r)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  r = Ray_2<R>(Point_2<R>(x1,y1),
                    Point_2<R>(x2,y2));
  w << r;
  return w;
}
#endif //CGAL_RAY_2_H




#ifdef TRIANGULATION_ELEMENT_2_H
template < class R >
Window_stream& operator<<(Window_stream &os,
                               const Triangulation_element_2<R> *t)
{

  switch(t->finite_vertices()){
  case 3 :
    {
      os <<  t->triangle();
      break;
    }
  case 2 :
    {
      int fi = 0;
      if (! t->is_finite(1)){
        fi = 1;
      }else if(! t->is_finite(2)){
        fi = 2;
      }
      os << Segment_2<R>((*t)[fi-1], (*t)[fi+1]);
      os << Ray_2<R>((*t)[fi+2],
                          (*t)[fi+2] + ((*t)[fi] - ORIGIN));
      os << Ray_2<R>((*t)[fi+1],
                          (*t)[fi+1] + ((*t)[fi] - ORIGIN));
   break;
    }
  case 1:
    {
      int fi = 0;
      if (t->is_finite(1)){
        fi = 1;
      }else if(t->is_finite(2)){
        fi = 2;
      }

      os << Ray_2<R>((*t)[fi],
                          (*t)[fi]+ ((*t)[fi+1] - ORIGIN));
      os << Ray_2<R>((*t)[fi],
                          (*t)[fi]+ ((*t)[fi+2] - ORIGIN));
      break;
    }
  }
  return os;
}
#endif // TRIANGULATION_ELEMENT_2_H

#ifdef CGAL_TRIANGULATION_2_H

template < class R >
Window_stream &operator<<(Window_stream &os,
                               Triangulation_2<R> &T)
{
  Triangulation_element_2<R> *t;

  forall_triangulation_elements(t, T){
    os << t;
  }
  return os;
}
#endif //  CGAL_TRIANGULATION_2_H


#ifdef PARABOLA_2_H
template< class R >
Window_stream& operator<<(Window_stream &w,
                               const Parabola_2<R> &par)
{
  double width = w.xmax()-w.xmin();
  double height = w.ymax()-w.ymin();
  double diag_sq = width*width + height*height;

  double lambda=1.0;
  
  Point_2<R> p1 = par.base(),
                  p2;

  while((par(lambda)-p1).operator*(par(lambda)-p1) < diag_sq){
    lambda *= 2.0;
  }

  while((par(lambda)-p1).operator*(par(lambda)-p1) > diag_sq){
    lambda /= 2.0;
  }

  lambda *= 2.0;
  double delta = lambda/50.0;

  int i;
  for(lambda = 0.0, i=0; i<=50 ; lambda += delta,i++){
    p2 = par(lambda);
    w << Segment_2<R>(p1,p2);
    p1 = p2;
  }
  p1 = par.base();
  for(lambda = 0.0, i=0; i<=50 ; lambda -= delta,i++){
    p2 = par(lambda);
    w << Segment_2<R>(p1,p2);
    p1 = p2;
  }
  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w,
                               Parabola_2<R> &par)
{
  Line_2<R> l;
  Point_2<R> p;
  w >> l;
  w >> p;
  par = Parabola_2<R>(l,p);
  w << par;
  return w;
}

#endif // PARABOLA_2_H

#ifdef PARABOLA_ARC_2_H
template< class R >
Window_stream& operator<<(Window_stream &w,
                               const Parabola_arc_2<R> &arc)
{
  double lambda, lmin, lmax;

  if (arc.lambda_target() > arc.lambda_source()){
    lmin = arc.lambda_source();
    lmax = arc.lambda_target();
  } else {
    lmin = arc.lambda_target();
    lmax = arc.lambda_source();
  }

  double delta = abs(lmax - lmin)/100.0;

  Point_2<R> p1 = arc.source();
  for(lambda = lmin;  lambda <= lmax; lambda += delta)
    {
      Point_2<R> p2 = arc.supporting_parabola()(lambda);
      w << Segment_2<R>(p1,p2);
      p1 = p2;
    }
  return w;
}
#endif // PARABOLA_ARC_2_H


#ifdef CGAL_CIRCLE_2_H
template< class R >
Window_stream& operator<<(Window_stream &w,
                               const Circle_2<R> &c)
{
  double cx = to_double(c.center().x()),
         cy = to_double(c.center().y()),
         r = to_double(c.squared_radius());

  w.draw_circle(cx, cy , sqrt(r));
  return w;
}

template< class R >
Window_stream& operator>>(Window_stream &w,
                               Circle_2<R> &c)
{
  double cx, cy, x1, y1;
  w.read_mouse(cx,cy);
  w.read_mouse_circle(cx,cy, x1, y1);
  Point_2<R> center(cx, cy),
                  p(x1, y1);

  Vector_2<R> v = center - p;
  R::FT sr = v*v;

  w.draw_circle(cx, cy , sqrt(sr));
  c = Circle_2<R>(center, sr);
  return w;
}

#endif // CGAL_CIRCLE_2_H



#ifdef CGAL_BBOX_2_H
inline Window_stream& operator<<(Window_stream &w,
                                      const Bbox_2 &b)
{
  line_style style = w.set_line_style(dotted);
  double xmin = b.xmin(),
         ymin = b.ymin(),
         xmax = b.xmax(),
         ymax = b.ymax();


  w.draw_segment(xmin, ymin, xmax, ymin);
  w.draw_segment(xmax, ymin, xmax, ymax);
  w.draw_segment(xmax, ymax, xmin, ymax);
  w.draw_segment(xmin, ymax, xmin, ymin);

  w.set_line_style(style);
  return w;
}
#endif // CGAL_BBOX_2_H

#endif // if 0

#endif // CGAL_WINDOW_STREAM_H

