/* 

Copyright (c) 1997 The CGAL Consortium

This software and related documentation is part of the 
Computational Geometry Algorithms Library (CGAL).

Permission to use, copy, and distribute this software and its 
documentation is hereby granted free of charge, provided that 
(1) it is not a component of a commercial product, and 
(2) this notice appears in all copies of the software and
    related documentation. 

CGAL may be distributed by any means, provided that the original
files remain intact, and no charge is made other than for
reasonable distribution costs.

CGAL may not be distributed as a component of any commercial
product without a prior license agreement with the authors.

This software and documentation is provided "as-is" and without 
warranty of any kind. In no event shall the CGAL Consortium be
liable for any damage of any kind.

The CGAL Consortium consists of Utrecht University (The Netherlands), 
ETH Zurich (Switzerland), Free University of Berlin (Germany), 
INRIA Sophia-Antipolis (France), Max-Planck-Institute Saarbrucken
(Germany), RISC Linz (Austria), and Tel-Aviv University (Israel).

*/

 
// Source: Window_stream.h
// Author: Andreas Fabri

#ifndef CGAL_Window_STREAM_H
#define CGAL_Window_STREAM_H


#include <CGAL/IO/Color.h>

#include <LEDA/window.h>
//#include <LEDA/REDEFINE_NAMES.h>
//typedef window leda_window;
//#include <LEDA/UNDEFINE_NAMES.h>

#include <CGAL/IO/optimisation_Window_stream.h>  // line added by Wieger
 
class CGAL_Window_stream : public leda_window  {

public:
  CGAL_Window_stream(int xpix = 300,
                     int ypix = 300)
    : leda_window(xpix, ypix)
  {
    set_frame_label("CGAL_Window_stream");
  }

  void init(const CGAL_Bbox_2 &b = CGAL_Bbox_2(0, 0, 1, 1))
  {
    leda_window::init(b.xmin(), b.xmax(), b.ymin());
    set_fg_color(leda_black);
  }

  void init(double xmin, double xmax, double ymin)
  {
    leda_window::init(xmin, xmax, ymin);
    set_fg_color(leda_black);
  }

  int read_mouse(double& x, double& y)
  {
    _button = leda_window::read_mouse(x,y);
    return _button;
  }

  int last_button_pressed()
  {
    return _button;
  }
  // we keep the foreground color in a field, as LEDA's window
  // can't store it in the shared library version
private:
  int _button;

};
 


#endif // CGAL_Window_STREAM_H

//  Each of the following operators is individually 
//  protected against multiple inclusion.

 
#ifdef CGAL_POINT_2_H
#ifndef CGAL_Window_STREAM_POINT_2
#define CGAL_Window_STREAM_POINT_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w, const CGAL_Point_2<R> &p)
{
  w.draw_point(CGAL_to_double(p.x()), CGAL_to_double(p.y()));
  return w;
}


template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w, CGAL_Point_2<R> &p)
{
  double x, y;
  w.read_mouse(x,y);
  w.draw_point(x,y);
  p = CGAL_Point_2<R>(R::RT(x), R::RT(y));
  return w;
}
#endif // CGAL_Window_STREAM_POINT_2
#endif //  CGAL_POINT_2_H
 

 
#ifdef CGAL_LINE_2_H
#ifndef CGAL_Window_STREAM_LINE_2
#define CGAL_Window_STREAM_LINE_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w,const CGAL_Segment_2<R> &s)
{
  w.draw_segment(CGAL_to_double(s.source().x()),
                 CGAL_to_double(s.source().y()),
                 CGAL_to_double(s.target().x()),
                 CGAL_to_double(s.target().y()));
  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w,
                               CGAL_Segment_2<R> &s)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  w.draw_segment(x1,y1, x2, y2);
  s = CGAL_Segment_2<R>(CGAL_Point_2<R>(R::FT(x1), R::FT(y1)),
                        CGAL_Point_2<R>(R::FT(x2), R::FT(y2)));
  return w;
}


template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w, const CGAL_Line_2<R> &l)
{
  CGAL_Point_2<R> p1 = l.point(),
                  p2 = p1 + l.direction().vector();

  w.draw_line(CGAL_to_double(p1.x()), CGAL_to_double(p1.y()),
              CGAL_to_double(p2.x()), CGAL_to_double(p2.y()));
  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w, CGAL_Line_2<R> &l)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  w.draw_line(x1,y1, x2, y2);

  l = CGAL_Line_2<R>(CGAL_Point_2<R>(x1,y1),
                     CGAL_Point_2<R>(x2,y2));
  return w;
}
#endif // CGAL_Window_STREAM_LINE_2
#endif //CGAL_LINE_2_H



#ifdef CGAL_RAY_2_H
#ifndef CGAL_Window_STREAM_RAY_2
#define CGAL_Window_STREAM_RAY_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w, const CGAL_Ray_2<R> &r)
{
  CGAL_Point_2<R> p = r.point(0);
  CGAL_Point_2<R> q = r.point(1);

  w.draw_ray(CGAL_to_double(p.x()),
             CGAL_to_double(p.y()),
             CGAL_to_double(q.x()),
             CGAL_to_double(q.y()));

  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w, CGAL_Ray_2<R> &r)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  r = CGAL_Ray_2<R>(CGAL_Point_2<R>(x1,y1),
                    CGAL_Point_2<R>(x2,y2));
  w << r;
  return w;
}
#endif // CGAL_Window_STREAM_RAY_2
#endif //CGAL_RAY_2_H

 


 
#ifdef CGAL_TRIANGLE_2_H
#ifndef CGAL_Window_STREAM_TRIANGLE_2
#define CGAL_Window_STREAM_TRIANGLE_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                               const CGAL_Triangle_2<R> &t)
{
  double x0 = CGAL_to_double(t.vertex(0).x()),
         y0 = CGAL_to_double(t.vertex(0).y()),
         x1 = CGAL_to_double(t.vertex(1).x()),
         y1 = CGAL_to_double(t.vertex(1).y()),
         x2 = CGAL_to_double(t.vertex(2).x()),
         y2 = CGAL_to_double(t.vertex(2).y());


  w.draw_segment(x0, y0, x1, y1);
  w.draw_segment(x1, y1, x2, y2);
  w.draw_segment(x2, y2, x0, y0);

  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w,
                               CGAL_Triangle_2<R> &t)
{
  double x0, y0, x1, y1, x2, y2;
  w.read_mouse(x0,y0);
  w.read_mouse_seg(x0,y0, x1, y1);
  w.draw_seg(x0,y0, x1, y1);
  w.read_mouse_seg(x1,y1, x2, y2);
  w.draw_seg(x1,y1, x2, y2);
  w.draw_seg(x2,y2, x0, y0);
  r = CGAL_Triangle_2<R>(CGAL_Point_2<R>(x0,y0),
                         CGAL_Point_2<R>(x1,y1),
                         CGAL_Point_2<R>(x2,y2));
  return w;
}
#endif // CGAL_Window_STREAM_TRIANGLE_2
#endif // CGAL_TRIANGLE_2_H
 


 
#ifdef CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_Window_STREAM_ISO_RECTANGLE_2
#define CGAL_Window_STREAM_ISO_RECTANGLE_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                               const CGAL_Iso_rectangle_2<R> &r)
{
  double xmin = CGAL_to_double(r.min().x()),
         ymin = CGAL_to_double(r.min().y()),
         xmax = CGAL_to_double(r.max().x()),
         ymax = CGAL_to_double(r.max().y());


  w.draw_segment(xmin, ymin, xmax, ymin);
  w.draw_segment(xmax, ymin, xmax, ymax);
  w.draw_segment(xmax, ymax, xmin, ymax);
  w.draw_segment(xmin, ymax, xmin, ymin);

  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w,
                               CGAL_Iso_rectangle_2<R> &r)
{
  double x1, y1, x2, y2;
  w.read_mouse(x1,y1);
  w.read_mouse_rect(x1,y1, x2, y2);
  r = CGAL_Iso_rectangle_2<R>(CGAL_Point_2<R>(x1,y1),
                              CGAL_Point_2<R>(x2,y2));
  w << r;
  return w;
}
#endif // CGAL_Window_STREAM_ISO_RECTANGLE_2
#endif // CGAL_ISO_RECTANGLE_2_H
 

 
#ifdef CGAL_PARABOLA_2_H
#ifndef CGAL_Window_STREAM_PARABOLA_2
#define CGAL_Window_STREAM_PARABOLA_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                               const CGAL_Parabola_2<R> &par)
{
  double width = w.xmax()-w.xmin();
  double height = w.ymax()-w.ymin();
  double diag_sq = width*width + height*height;

  double lambda=1.0;
  
  CGAL_Point_2<R> p1 = par.base(),
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
    w << CGAL_Segment_2<R>(p1,p2);
    p1 = p2;
  }
  p1 = par.base();
  for(lambda = 0.0, i=0; i<=50 ; lambda -= delta,i++){
    p2 = par(lambda);
    w << CGAL_Segment_2<R>(p1,p2);
    p1 = p2;
  }
  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w,
                               CGAL_Parabola_2<R> &par)
{
  CGAL_Line_2<R> l;
  CGAL_Point_2<R> p;
  w >> l;
  w >> p;
  par = CGAL_Parabola_2<R>(l,p);
  w << par;
  return w;
}
#endif // CGAL_Window_STREAM_PARABOLA_2
#endif // CGAL_PARABOLA_2_H
 

 
#ifdef CGAL_PARABOLA_ARC_2_H
#ifndef CGAL_Window_STREAM_PARABOLA_ARC_2
#define CGAL_Window_STREAM_PARABOLA_ARC_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                               const CGAL_Parabola_arc_2<R> &arc)
{
  double lambda, lmin, lmax;

  if (arc.lambda_target() > arc.lambda_source()){
    lmin = arc.lambda_source();
    lmax = arc.lambda_target();
  } else {
    lmin = arc.lambda_target();
    lmax = arc.lambda_source();
  }

  double delta = CGAL_abs(lmax - lmin)/100.0;

  CGAL_Point_2<R> p1 = arc.source();
  for(lambda = lmin;  lambda <= lmax; lambda += delta)
    {
      CGAL_Point_2<R> p2 = arc.supporting_parabola()(lambda);
      w << CGAL_Segment_2<R>(p1,p2);
      p1 = p2;
    }
  return w;
}
#endif //CGAL_Window_STREAM_PARABOLA_ARC_2
#endif // CGAL_PARABOLA_ARC_2_H
 


 
#ifdef CGAL_CIRCLE_2_H
#ifndef CGAL_Window_STREAM_CIRCLE_2
#define CGAL_Window_STREAM_CIRCLE_2
template< class R >
CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                               const CGAL_Circle_2<R> &c)
{
  double cx = CGAL_to_double(c.center().x()),
         cy = CGAL_to_double(c.center().y()),
         r = CGAL_to_double(c.squared_radius());

  w.draw_circle(cx, cy , sqrt(r));
  return w;
}

template< class R >
CGAL_Window_stream& operator>>(CGAL_Window_stream &w,
                               CGAL_Circle_2<R> &c)
{
  double cx, cy, x1, y1;
  w.read_mouse(cx,cy);
  w.read_mouse_circle(cx,cy, x1, y1);
  CGAL_Point_2<R> center(cx, cy),
                  p(x1, y1);

  CGAL_Vector_2<R> v = center - p;
  R::FT sr = v*v;

  w.draw_circle(cx, cy , sqrt(sr));
  c = CGAL_Circle_2<R>(center, sr);
  return w;
}
#endif // CGAL_Window_STREAM_CIRCLE_2
#endif // CGAL_CIRCLE_2_H
 


 
#ifndef CGAL_Window_STREAM_COLOR_2
#define CGAL_Window_STREAM_COLOR_2
inline CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                                      const CGAL_Color& c)
{
  w.set_fg_color(leda_color(c.red(), c.green(), c.blue()));

  return w;
}
#endif // CGAL_Window_STREAM_COLOR_2
 

 
#ifdef CGAL_BBOX_2_H
#ifndef CGAL_Window_STREAM_BBOX_2
#define CGAL_Window_STREAM_BBOX_2
inline CGAL_Window_stream& operator<<(CGAL_Window_stream &w,
                                      const CGAL_Bbox_2 &b)
{
  line_style style = w.set_line_style(leda_dotted);
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
#endif // CGAL_Window_STREAM_BBOX_2
#endif // CGAL_BBOX_2_H
 


 
#ifdef CGAL_TRIANGULATION_2_H
#ifndef CGAL_Window_STREAM_TRIANGULATION_2_H
#define CGAL_Window_STREAM_TRIANGULATION_2_H
template < class I >
CGAL_Window_stream&
operator<<(CGAL_Window_stream& os,
           const CGAL_Triangulation_2<I> &T)
{
   CGAL_Triangulation_2<I>::Edge_iterator it = T.edges_begin();

    while(it != T.edges_end()){
        os << T.segment(it);
        ++it;
    }

    return os;
}
#endif // CGAL_Window_STREAM_TRIANGULATION_2_H
#endif // CGAL_TRIANGULATION_2_H


#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H
#ifndef CGAL_Window_STREAM_DELAUNAY_TRIANGULATION_2_H
#define CGAL_Window_STREAM_DELAUNAY_TRIANGULATION_2_H
template < class I >
CGAL_Window_stream&
operator<<(CGAL_Window_stream& os,
           const CGAL_Delaunay_triangulation_2<I> &T)
{
   CGAL_Delaunay_triangulation_2<I>::Edge_iterator it = T.edges_begin();

    while(it != T.edges_end()){
        os << T.segment(it);
        ++it;
    }

    return os;
}
#endif // CGAL_Window_STREAM_DELAUNAY_TRIANGULATION_2_H
#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

