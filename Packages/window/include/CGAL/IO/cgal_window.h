// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
// 
// file          : include/CGAL/IO/cgal_window.h
// package       : window (2.8.5)
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// revision      : 2.8.5
// revision_date : 22 June 2001 
// author(s)     : Matthias Baesken
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_LEDA_WINDOW_H
#define CGAL_LEDA_WINDOW_H

#include <CGAL/IO/Color.h>
#include <CGAL/LEDA/window.h>
#include <CGAL/IO/esprit_logo.xpm>


CGAL_BEGIN_NAMESPACE


typedef CGAL::window        Window_stream;

inline
CGAL::window&
operator<<(CGAL::window& w, const Color& c)
{
  w.set_fg_color(CGAL::color(c.red(), c.green(), c.blue()));
  return w;
}


inline
void
cgalize(CGAL::window& w)
{
  w.set_frame_label("CGAL-2.3");
  w.set_icon_label("CGAL");
  w.set_line_width( 2);
  w.set_icon_pixrect( w.create_pixrect((const char**) esprit_logo));
}

inline
CGAL::window*
create_demo_window( float w = 512.0, float h = 512.0,
                         const char* str = "CGAL",
                         double x_extension = 1.0)
{
  CGAL::window* Wptr = new CGAL::window((int) w, (int) h);
  cgalize( *Wptr);
  double y_extension = x_extension * h / w;
  Wptr->init(-x_extension, x_extension, -y_extension);
  Wptr->set_frame_label( str);
  return Wptr;
}


inline
CGAL::window*
create_and_display_demo_window(float w = 512.0, float h = 512.0,
                                    const char* str = "CGAL",
                                    double x_extension = 1.0)
{
  CGAL::window* Wptr = new CGAL::window((int) w, (int) h);
  cgalize( *Wptr);
  double y_extension = x_extension * h / w;
  Wptr->init(-x_extension, x_extension, -y_extension);
  Wptr->set_frame_label( str);
  Wptr->display();
  return Wptr;
}

CGAL_END_NAMESPACE


#endif // CGAL_LEDA_WINDOW_H

CGAL_BEGIN_NAMESPACE

//  Each of the following operators is individually
//  protected against multiple inclusion.

#ifdef CGAL_POINT_2_H
#ifndef CGAL_LEDA_WINDOW_POINT_2
#define CGAL_LEDA_WINDOW_POINT_2
template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Point_2<R>& p)
{
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  w.draw_point(x,y);
  
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Point_2<R>& p)
{
  typedef typename R::RT RT;
  CGAL::window_point l_p;
  drawing_mode save = w.set_mode(xor_mode);
  if (w.read(l_p))
  {
      double x = l_p.xcoord();
      double y = l_p.ycoord();
      w.draw_point(x,y);
      w.set_mode( save);
      w.draw_point(x,y);
      
      p = Point_2<R>( RT(x), RT(y));
  }
  else
  {
      w.set_mode( save);
  }
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Point_2<R>& p)
{
  typedef typename R::RT RT;
  CGAL::window_point l_p;
  if (w.read(l_p))
  {
      double x = l_p.xcoord();
      double y = l_p.ycoord();
      p = Point_2<R>( RT(x), RT(y));
  }
  return w;
}

template <class R>
void
read_mouse_plus(CGAL::window& w, Point_2<R>& p, int& button)
{
  typedef typename R::RT RT;
  double x, y;
  button = w.read_mouse(x,y);

  w.draw_point(x,y);
  
  p = Point_2<R>(RT(x), RT(y));
}

#endif // CGAL_LEDA_WINDOW_POINT_2
#endif // CGAL_POINT_2_H


#ifdef CGAL_SEGMENT_2_H
#ifndef CGAL_LEDA_WINDOW_SEGMENT_2
#define CGAL_LEDA_WINDOW_SEGMENT_2
template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Segment_2<R>& s)
{
  w.draw_segment(CGAL::to_double(s.source().x()),
                 CGAL::to_double(s.source().y()),
                 CGAL::to_double(s.target().x()),
                 CGAL::to_double(s.target().y()));
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Segment_2<R>& s)
{
  typedef  typename R::RT  RT;
  CGAL::window_point l_p1, l_p2;
  drawing_mode save = w.set_mode(xor_mode);
  if ( w.read( l_p1, l_p2))
  {
      double x1 = l_p1.xcoord();
      double y1 = l_p1.ycoord();
      double x2 = l_p2.xcoord();
      double y2 = l_p2.ycoord();
      w.set_mode( save);
      w.draw_segment(x1,y1, x2, y2);
      s = Segment_2<R>(Point_2<R>( RT(x1), RT(y1)),
                            Point_2<R>( RT(x2), RT(y2)));
  }
  else
  {
      w.set_mode( save);
  }
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Segment_2<R>& s)
{
  typedef  typename R::RT  RT;
  CGAL::window_point l_p1, l_p2;
  
  if ( w.read( l_p1, l_p2))
  {
      double x1 = l_p1.xcoord();
      double y1 = l_p1.ycoord();
      double x2 = l_p2.xcoord();
      double y2 = l_p2.ycoord();
      s = Segment_2<R>(Point_2<R>( RT(x1), RT(y1)),
                            Point_2<R>( RT(x2), RT(y2)));
  }
  return w;
}
#endif // CGAL_LEDA_WINDOW_SEGMENT_2
#endif // CGAL_SEGMENT_2_H


#ifdef CGAL_LINE_2_H
#ifndef CGAL_LEDA_WINDOW_LINE_2
#define CGAL_LEDA_WINDOW_LINE_2

template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Line_2<R>& l)
{
  Point_2<R> p1 = l.point(),
                  p2 = p1 + l.direction().vector();
  w.draw_line(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),
              CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Line_2<R>& l)
{
  typedef  typename R::RT  RT;
  CGAL::window_point l_p1, l_p2;

  drawing_mode save = w.set_mode(xor_mode);
  if ( w.read( l_p1, l_p2))
  {
      double x1 = l_p1.xcoord();
      double y1 = l_p1.ycoord();
      double x2 = l_p2.xcoord();
      double y2 = l_p2.ycoord();
      w.set_mode( save);
      w.draw_line(x1,y1, x2, y2);
      l = Line_2<R>(Point_2<R>( RT(x1), RT(y1)),
                         Point_2<R>( RT(x2), RT(y2)));
  }
  else
  {
      w.set_mode( save);
  }
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Line_2<R>& l)
{
  typedef  typename R::RT  RT;
  CGAL::window_point l_p1, l_p2;  
  if ( w.read( l_p1,l_p2))
  {
      double x1 = l_p1.xcoord();
      double y1 = l_p1.ycoord();
      double x2 = l_p2.xcoord();
      double y2 = l_p2.ycoord();
      l = Line_2<R>(Point_2<R>( RT(x1), RT(y1)),
                         Point_2<R>( RT(x2), RT(y2)));
  }
  return w;
}
#endif // CGAL_LEDA_WINDOW_LINE_2
#endif //CGAL_LINE_2_H

#ifdef CGAL_RAY_2_H
#ifndef CGAL_LEDA_WINDOW_RAY_2
#define CGAL_LEDA_WINDOW_RAY_2
template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Ray_2<R>& r)
{
  Point_2<R> p = r.point(0);
  Point_2<R> q = r.point(1);
  w.draw_ray(CGAL::to_double(p.x()),
             CGAL::to_double(p.y()),
             CGAL::to_double(q.x()),
             CGAL::to_double(q.y()));
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Ray_2<R>& r)
{
  typedef  typename R::RT  RT;
  CGAL::window_point l_p1, l_p2;
  drawing_mode save = w.set_mode(xor_mode);
  if ( w.read( l_p1, l_p2))
  {
      double x1 = l_p1.xcoord();
      double y1 = l_p1.ycoord();
      double x2 = l_p2.xcoord();
      double y2 = l_p2.ycoord();
      w.set_mode( save);
      r = Ray_2<R>(Point_2<R>( RT(x1), RT(y1)),
                        Point_2<R>( RT(x2), RT(y2)));
      w << r;
  }
  else
  {
      w.set_mode( save);
  }
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Ray_2<R>& r)
{
  typedef  typename R::RT  RT;
  CGAL::window_point l_p1, l_p2;
  if ( w.read( l_p1, l_p2))
  {
      double x1 = l_p1.xcoord();
      double y1 = l_p1.ycoord();
      double x2 = l_p2.xcoord();
      double y2 = l_p2.ycoord();
      r = Ray_2<R>(Point_2<R>( RT(x1), RT(y1)),
                        Point_2<R>( RT(x2), RT(y2)));
  }
  return w;
}
#endif // CGAL_LEDA_WINDOW_RAY_2
#endif //CGAL_RAY_2_H

#ifdef CGAL_TRIANGLE_2_H
#ifndef CGAL_LEDA_WINDOW_TRIANGLE_2
#define CGAL_LEDA_WINDOW_TRIANGLE_2
template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Triangle_2<R>& t)
{
  double x0 = CGAL::to_double(t.vertex(0).x()),
         y0 = CGAL::to_double(t.vertex(0).y()),
         x1 = CGAL::to_double(t.vertex(1).x()),
         y1 = CGAL::to_double(t.vertex(1).y()),
         x2 = CGAL::to_double(t.vertex(2).x()),
         y2 = CGAL::to_double(t.vertex(2).y());
	 
  CGAL::color cl = w.get_fill_color();
  
  if (cl != CGAL::invisible) { // draw filled triangle ...
    w.draw_filled_triangle(CGAL::window_point(x0,y0), CGAL::window_point(x1,y1), CGAL::window_point(x2,y2), cl); 
  }
	 
  w.draw_segment(x0, y0, x1, y1);
  w.draw_segment(x1, y1, x2, y2);
  w.draw_segment(x2, y2, x0, y0);
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Triangle_2<R>& t)
{
  typedef typename R::RT   RT;
  double x,y;
  int key = 0;

  w.set_state( 1);

  CGAL::window_point first, second, third;
  drawing_mode save = w.set_mode(xor_mode);
  
  if ( !( w.read(first))) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_seg( first.x(), first.y(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  else
  {
      w.draw_segment( first.x(), first.y(), x, y);
      second = CGAL::window_point( x, y);
  }
  key = w.read_mouse_seg( second.x(), second.y(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.draw_segment( first.x(), first.y(), second.x(), second.y() );
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  else
  {
      w.draw_segment( second.x(), second.y(), x, y);
      third = CGAL::window_point( x, y);
  }
  w.draw_segment( first.x(), first.y(), second.x(), second.y() );
  w.draw_segment( second.x(), second.y(), third.x(), third.y() );
  double x0 = first.x();
  double y0 = first.y();
  double x1 = second.x();
  double y1 = second.y();
  double x2 = third.x();
  double y2 = third.y();
  w.set_mode( save);
  
  CGAL::color cl = w.get_fill_color();
  
  if (cl != CGAL::invisible) { // draw filled triangle ...
    w.draw_filled_triangle(CGAL::window_point(x0,y0), CGAL::window_point(x1,y1), CGAL::window_point(x2,y2), cl); 
  }  
  
  w.draw_segment(x0,y0, x1, y1);
  w.draw_segment(x1,y1, x2, y2);
  w.draw_segment(x2,y2, x0, y0);

  t = Triangle_2<R>(Point_2<R>( RT(x0), RT(y0)),
                         Point_2<R>( RT(x1), RT(y1)),
                         Point_2<R>( RT(x2), RT(y2)));
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Triangle_2<R>& t)
{
  typedef typename R::RT   RT;
  double x,y;
  int key = 0;

  w.set_state( 1);

  CGAL::window_point first, second, third;
  drawing_mode save = w.set_mode(xor_mode);
  if ( !( w.read(first))) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_seg( first.x(), first.y(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  else
  {
      w.draw_segment( first.x(), first.y(), x, y);
      second = CGAL::window_point( x, y);
  }
  key = w.read_mouse_seg( second.x(), second.y(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.draw_segment( first.x(), first.y(), second.x(), second.y() );
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  else
  {
      w.draw_segment( second.x(), second.y(), x, y);
      third = CGAL::window_point( x, y);
  }
  w.draw_segment( first.x(), first.y(), second.x(), second.y() );
  w.draw_segment( second.x(), second.y(), third.x(), third.y() );
  double x0 = first.x();
  double y0 = first.y();
  double x1 = second.x();
  double y1 = second.y();
  double x2 = third.x();
  double y2 = third.y();
  w.set_mode( save);
  t = Triangle_2<R>(Point_2<R>( RT(x0), RT(y0)),
                         Point_2<R>( RT(x1), RT(y1)),
                         Point_2<R>( RT(x2), RT(y2)));
  return w;
}


#endif // CGAL_LEDA_WINDOW_TRIANGLE_2
#endif // CGAL_TRIANGLE_2_H

#ifdef CGAL_CIRCLE_2_H
#ifndef CGAL_LEDA_WINDOW_CIRCLE_2
#define CGAL_LEDA_WINDOW_CIRCLE_2
template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Circle_2<R>& c)
{
  double cx = CGAL::to_double(c.center().x()),
         cy = CGAL::to_double(c.center().y()),
         r  = CGAL::to_double(c.squared_radius());
	 
  CGAL::color cl = w.get_fill_color();
  
  if (cl != CGAL::invisible) { // draw filled circle ...
    w.draw_disc(cx, cy , std::sqrt(r), cl); 
  }	 
	 
  w.draw_circle(cx, cy , std::sqrt(r));
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Circle_2<R>& c)
{
  typedef  typename R::RT  RT;
  double x,y;
  int key = 0;

  w.set_state( 1);

  CGAL::window_point p;
  drawing_mode save = w.set_mode(xor_mode);
  if ( !( w.read( p))) { w.set_mode( save); return w; }
  
  double cx = p.x();
  double cy = p.y();
  Point_2<R> center = Point_2<R>( RT(cx), RT(cy));
  
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_circle(cx, cy, x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);
      return w;
  }
  double dx = x - cx;
  double dy = y - cy;
  double sqr = dx*dx+dy*dy;
  w.set_mode( save);
  w.set_buttons( save_but);

  CGAL::color cl = w.get_fill_color();
  
  if (cl != CGAL::invisible) { // draw filled circle ...
    w.draw_disc(cx, cy , std::sqrt(sqr), cl); 
  }  
  
  w.draw_circle(cx, cy , std::sqrt(sqr));
  c = Circle_2<R>(center, RT(sqr));
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Circle_2<R>& c)
{
  typedef  typename R::RT  RT;
  double x,y;
  int key = 0;

  w.set_state( 1);

  CGAL::window_point p;
  drawing_mode save = w.set_mode(xor_mode);
  if ( !( w.read( p))) { w.set_mode( save); return w; }
  double cx = p.x();
  double cy = p.y();
  Point_2<R> center = Point_2<R>( RT(cx), RT(cy));
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_circle(cx, cy, x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  double dx = x - cx;
  double dy = y - cy;
  double sqr = dx*dx+dy*dy;
  w.set_mode( save);
  w.set_buttons( save_but);
  c = Circle_2<R>(center, RT(sqr));
  return w;
}
#endif // CGAL_LEDA_WINDOW_CIRCLE_2
#endif // CGAL_CIRCLE_2_H

#ifdef CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_LEDA_WINDOW_ISO_RECTANGLE_2
#define CGAL_LEDA_WINDOW_ISO_RECTANGLE_2
template< class R >
CGAL::window&
operator<<(CGAL::window& w, const Iso_rectangle_2<R>& r)
{
  double xmin = CGAL::to_double(r.min().x()),
         ymin = CGAL::to_double(r.min().y()),
         xmax = CGAL::to_double(r.max().x()),
         ymax = CGAL::to_double(r.max().y());

  CGAL::color cl = w.get_fill_color();
  
  if (cl != CGAL::invisible) { // draw filled rectangle ...
    w.draw_filled_rectangle(xmin, ymin, xmax, ymax, cl); 
  }		 
	 
  w.draw_segment(xmin, ymin, xmax, ymin);
  w.draw_segment(xmax, ymin, xmax, ymax);
  w.draw_segment(xmax, ymax, xmin, ymax);
  w.draw_segment(xmin, ymax, xmin, ymin);
  return w;
}

template< class R >
CGAL::window&
operator>>(CGAL::window& w, Iso_rectangle_2<R>& r)
{
  typedef typename R::RT    RT;
  double x,y;
  int key = 0;

  w.set_state( 1);

  CGAL::window_point first;
  drawing_mode save = w.set_mode(xor_mode);
  
  if ( !( w.read( first))) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_rect( first.x(), first.y(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  r = Iso_rectangle_2<R>(Point_2<R>( RT(first.x()),
                                               RT(first.y())),
                              Point_2<R>( RT(x), RT(y)));
  w.set_mode( save);

  CGAL::color cl = w.get_fill_color();
  
  if (cl != CGAL::invisible) { // draw filled rectangle ...
    w.draw_filled_rectangle(first.x(), first.y(), x, y, cl); 
  }	  
  
  w.draw_rectangle( first.x(), first.y(), x, y);
  w.set_buttons( save_but);
  return w;
}

template< class R >
CGAL::window&
read(CGAL::window& w, Iso_rectangle_2<R>& r)
{
  typedef typename R::RT    RT;
  double x,y;
  int key = 0;

  w.set_state( 1);

  CGAL::window_point first;
  drawing_mode save = w.set_mode(xor_mode);
  if ( !( w.read( first))) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_rect( first.x(), first.y(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);

      w.set_state( 0);

      return w;
  }
  r = Iso_rectangle_2<R>(Point_2<R>( RT(first.x()),
                                               RT(first.y())),
                              Point_2<R>( RT(x), RT(y)));
  w.set_mode( save);
  w.set_buttons( save_but);
  return w;
}
#endif // CGAL_LEDA_WINDOW_ISO_RECTANGLE_2
#endif // CGAL_ISO_RECTANGLE_2_H

#ifdef CGAL_BBOX_2_H
#ifndef CGAL_LEDA_WINDOW_BBOX_2
#define CGAL_LEDA_WINDOW_BBOX_2
inline
CGAL::window&
operator<<(CGAL::window& w, const Bbox_2& b)
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
#endif // CGAL_LEDA_WINDOW_BBOX_2
#endif // CGAL_BBOX_2_H


CGAL_END_NAMESPACE

#ifndef IO_TRIANGULATION_WINDOW_STREAM_H
#include <CGAL/IO/triangulation_Window_stream.h>
#endif  // IO_TRIANGULATION_WINDOW_STREAM_H
#ifndef IO_OPTIMISATION_WINDOW_STREAM_H
#include <CGAL/IO/optimisation_Window_stream.h>
#endif // IO_OPTIMISATION_WINDOW_STREAM_H
#ifndef IO_POLYGON_WINDOW_STREAM_H
#include <CGAL/IO/polygon_Window_stream.h>
#endif // IO_POLYGON_WINDOW_STREAM_H

