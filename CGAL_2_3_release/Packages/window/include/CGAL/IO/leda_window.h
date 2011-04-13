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
// release       : $CGAL_Revision: CGAL-2.3-I-73 $
// release_date  : $CGAL_Date: 2001/06/19 $
// 
// file          : include/CGAL/IO/leda_window.h
// package       : window (2.8.4)
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
// revision      : 2.8.4
// revision_date : 21 June 2001
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_LEDA_WINDOW_H
#define CGAL_LEDA_WINDOW_H

#include <CGAL/IO/Color.h>
#include <LEDA/window.h>
// #include <CGAL/IO/cgal_icon.xpm>
#include <CGAL/IO/esprit_logo.xpm>


CGAL_BEGIN_NAMESPACE


// #if ( __LEDA__ < 360 )
// #define leda_window window
// #define leda_color  color
// #endif

#if ( __LEDA__ < 361 )
#define leda_black black
#define leda_dotted dotted
#define leda_xor_mode xor_mode
#define leda_blue blue
#define leda_red red
#endif

#if (( __LEDA__ < 362 ) || defined(NOTEBOOK))
#define  leda_drawing_mode  drawing_mode
#endif

typedef leda_window        Window_stream;

inline
leda_window&
operator<<(leda_window& w, const Color& c)
{
  w.set_fg_color(leda_color(c.red(), c.green(), c.blue()));
  return w;
}


inline
void
cgalize(leda_window& w)
{
  w.set_frame_label("CGAL-2.3");
  w.set_icon_label("CGAL");
  w.set_line_width( 2);
  w.set_icon_pixrect( w.create_pixrect( esprit_logo));
}

inline
leda_window*
create_demo_window( float w = 512.0, float h = 512.0,
                         const char* str = "CGAL",
                         double x_extension = 1.0)
{
  leda_window* Wptr = new leda_window( w, h);
  cgalize( *Wptr);
  double y_extension = x_extension * h / w;
  Wptr->init(-x_extension, x_extension, -y_extension);
  Wptr->set_frame_label( str);
  return Wptr;
}


inline
leda_window*
create_and_display_demo_window(float w = 512.0, float h = 512.0,
                                    const char* str = "CGAL",
                                    double x_extension = 1.0)
{
  leda_window* Wptr = new leda_window( w, h);
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
leda_window&
operator<<(leda_window& w, const Point_2<R>& p)
{
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  w.draw_point(x,y);
  
  return w;
}

template< class R >
leda_window&
operator>>(leda_window& w, Point_2<R>& p)
{
  typedef typename R::RT RT;
  leda_point l_p;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if (w >> l_p)
  {
      double x = l_p.xcoord();
      double y = l_p.ycoord();
      w << l_p;
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
leda_window&
read(leda_window& w, Point_2<R>& p)
{
  typedef typename R::RT RT;
  leda_point l_p;
  if (w >> l_p)
  {
      double x = l_p.xcoord();
      double y = l_p.ycoord();
      p = Point_2<R>( RT(x), RT(y));
  }
  return w;
}

template <class R>
void
read_mouse_plus(leda_window& w, Point_2<R>& p, int& button)
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
leda_window&
operator<<(leda_window& w, const Segment_2<R>& s)
{
  w.draw_segment(CGAL::to_double(s.source().x()),
                 CGAL::to_double(s.source().y()),
                 CGAL::to_double(s.target().x()),
                 CGAL::to_double(s.target().y()));
  return w;
}

template< class R >
leda_window&
operator>>(leda_window& w, Segment_2<R>& s)
{
  typedef  typename R::RT  RT;
  leda_segment l_s;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( w.read( l_s))
  {
      double x1 = l_s.xcoord1();
      double y1 = l_s.ycoord1();
      double x2 = l_s.xcoord2();
      double y2 = l_s.ycoord2();
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
leda_window&
read(leda_window& w, Segment_2<R>& s)
{
  typedef  typename R::RT  RT;
  leda_segment l_s;
  if ( w.read( l_s))
  {
      double x1 = l_s.xcoord1();
      double y1 = l_s.ycoord1();
      double x2 = l_s.xcoord2();
      double y2 = l_s.ycoord2();
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
leda_window&
operator<<(leda_window& w, const Line_2<R>& l)
{
  Point_2<R> p1 = l.point(),
                  p2 = p1 + l.direction().vector();
  w.draw_line(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),
              CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));
  return w;
}

template< class R >
leda_window&
operator>>(leda_window& w, Line_2<R>& l)
{
  typedef  typename R::RT  RT;
  leda_segment l_s;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( w.read( l_s))
  {
      double x1 = l_s.xcoord1();
      double y1 = l_s.ycoord1();
      double x2 = l_s.xcoord2();
      double y2 = l_s.ycoord2();
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
leda_window&
read(leda_window& w, Line_2<R>& l)
{
  typedef  typename R::RT  RT;
  leda_segment l_s;
  if ( w.read( l_s))
  {
      double x1 = l_s.xcoord1();
      double y1 = l_s.ycoord1();
      double x2 = l_s.xcoord2();
      double y2 = l_s.ycoord2();
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
leda_window&
operator<<(leda_window& w, const Ray_2<R>& r)
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
leda_window&
operator>>(leda_window& w, Ray_2<R>& r)
{
  typedef  typename R::RT  RT;
  leda_segment l_s;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( w.read( l_s))
  {
      double x1 = l_s.xcoord1();
      double y1 = l_s.ycoord1();
      double x2 = l_s.xcoord2();
      double y2 = l_s.ycoord2();
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
leda_window&
read(leda_window& w, Ray_2<R>& r)
{
  typedef  typename R::RT  RT;
  leda_segment l_s;
  if ( w.read( l_s))
  {
      double x1 = l_s.xcoord1();
      double y1 = l_s.ycoord1();
      double x2 = l_s.xcoord2();
      double y2 = l_s.ycoord2();
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
leda_window&
operator<<(leda_window& w, const Triangle_2<R>& t)
{
  double x0 = CGAL::to_double(t.vertex(0).x()),
         y0 = CGAL::to_double(t.vertex(0).y()),
         x1 = CGAL::to_double(t.vertex(1).x()),
         y1 = CGAL::to_double(t.vertex(1).y()),
         x2 = CGAL::to_double(t.vertex(2).x()),
         y2 = CGAL::to_double(t.vertex(2).y());

  leda_color cl = w.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled triangle ...
    w.draw_filled_triangle(leda_point(x0,y0), leda_point(x1,y1), leda_point(x2,y2), cl); 
  }	 
	 
	 
  w.draw_segment(x0, y0, x1, y1);
  w.draw_segment(x1, y1, x2, y2);
  w.draw_segment(x2, y2, x0, y0);
  return w;
}

template< class R >
leda_window&
operator>>(leda_window& w, Triangle_2<R>& t)
{
  typedef typename R::RT   RT;
  double x,y;
  int key = 0;
#if ( __LEDA__ < 362 )
  w.state = 1;
#else
  w.set_state( 1);
#endif // __LEDA__ < ...
  leda_point first, second, third, p;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( !( w >> first)) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_seg( first.xcoord(), first.ycoord(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  else
  {
      w << leda_segment( first.xcoord(), first.ycoord(), x, y);
      second = leda_point( x, y);
  }
  key = w.read_mouse_seg( second.xcoord(), second.ycoord(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w << leda_segment( first, second );
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  else
  {
      w << leda_segment( second.xcoord(), second.ycoord(), x, y);
      third = leda_point( x, y);
  }
  w << leda_segment( first, second );
  w << leda_segment( second, third );
  double x0 = first.xcoord();
  double y0 = first.ycoord();
  double x1 = second.xcoord();
  double y1 = second.ycoord();
  double x2 = third.xcoord();
  double y2 = third.ycoord();
  w.set_mode( save);
  
  leda_color cl = w.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled triangle ...
    w.draw_filled_triangle(leda_point(x0,y0), leda_point(x1,y1), leda_point(x2,y2), cl); 
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
leda_window&
read(leda_window& w, Triangle_2<R>& t)
{
  typedef typename R::RT   RT;
  double x,y;
  int key = 0;
#if ( __LEDA__ < 362 )
  w.state = 1;
#else
  w.set_state( 1);
#endif // __LEDA__ < ...
  leda_point first, second, third, p;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( !( w >> first)) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_seg( first.xcoord(), first.ycoord(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  else
  {
      w << leda_segment( first.xcoord(), first.ycoord(), x, y);
      second = leda_point( x, y);
  }
  key = w.read_mouse_seg( second.xcoord(), second.ycoord(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w << leda_segment( first, second );
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  else
  {
      w << leda_segment( second.xcoord(), second.ycoord(), x, y);
      third = leda_point( x, y);
  }
  w << leda_segment( first, second );
  w << leda_segment( second, third );
  double x0 = first.xcoord();
  double y0 = first.ycoord();
  double x1 = second.xcoord();
  double y1 = second.ycoord();
  double x2 = third.xcoord();
  double y2 = third.ycoord();
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
leda_window&
operator<<(leda_window& w, const Circle_2<R>& c)
{
  double cx = CGAL::to_double(c.center().x()),
         cy = CGAL::to_double(c.center().y()),
         r  = CGAL::to_double(c.squared_radius());
	 
  leda_color cl = w.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled circle ...
    w.draw_disc(cx, cy , std::sqrt(r), cl); 
  }		 
	 
  w.draw_circle(cx, cy , std::sqrt(r));
  return w;
}

template< class R >
leda_window&
operator>>(leda_window& w, Circle_2<R>& c)
{
  typedef  typename R::RT  RT;
  double x,y;
  int key = 0;
#if ( __LEDA__ < 362 )
  w.state = 1;
#else
  w.set_state( 1);
#endif // __LEDA__ < ...
  leda_point p;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( !( w.read( p))) { w.set_mode( save); return w; }
  double cx = p.xcoord();
  double cy = p.ycoord();
  Point_2<R> center = Point_2<R>( RT(cx), RT(cy));
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_circle(cx, cy, x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  double dx = x - cx;
  double dy = y - cy;
  double sqr = dx*dx+dy*dy;
  w.set_mode( save);
  w.set_buttons( save_but);
  
  leda_color cl = w.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled circle ...
    w.draw_disc(cx, cy , std::sqrt(sqr), cl); 
  }    
  
  w.draw_circle(cx, cy , std::sqrt(sqr));
  c = Circle_2<R>(center, RT(sqr));
  return w;
}

template< class R >
leda_window&
read(leda_window& w, Circle_2<R>& c)
{
  typedef  typename R::RT  RT;
  double x,y;
  int key = 0;
#if ( __LEDA__ < 362 )
  w.state = 1;
#else
  w.set_state( 1);
#endif // __LEDA__ < ...
  leda_point p;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( !( w.read( p))) { w.set_mode( save); return w; }
  double cx = p.xcoord();
  double cy = p.ycoord();
  Point_2<R> center = Point_2<R>( RT(cx), RT(cy));
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_circle(cx, cy, x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
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
leda_window&
operator<<(leda_window& w, const Iso_rectangle_2<R>& r)
{
  double xmin = CGAL::to_double(r.min().x()),
         ymin = CGAL::to_double(r.min().y()),
         xmax = CGAL::to_double(r.max().x()),
         ymax = CGAL::to_double(r.max().y());
	 
  leda_color cl = w.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled rectangle ...
    w.draw_filled_rectangle(xmin, ymin, xmax, ymax, cl); 
  }	 
	 
  w.draw_segment(xmin, ymin, xmax, ymin);
  w.draw_segment(xmax, ymin, xmax, ymax);
  w.draw_segment(xmax, ymax, xmin, ymax);
  w.draw_segment(xmin, ymax, xmin, ymin);
  return w;
}

template< class R >
leda_window&
operator>>(leda_window& w, Iso_rectangle_2<R>& r)
{
  typedef typename R::RT    RT;
  double x,y;
  int key = 0;
#if ( __LEDA__ < 362 )
  w.state = 1;
#else
  w.set_state( 1);
#endif // __LEDA__ < ...
  leda_point first, second;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( !( w.read( first))) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_rect( first.xcoord(), first.ycoord(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  r = Iso_rectangle_2<R>(Point_2<R>( RT(first.xcoord()),
                                               RT(first.ycoord())),
                              Point_2<R>( RT(x), RT(y)));
  w.set_mode( save);
  
  leda_color cl = w.get_fill_color();
  
  if (cl != leda_invisible) { // draw filled rectangle ...
    w.draw_filled_rectangle(first.xcoord(), first.ycoord(), x, y, cl); 
  }  
  
  w.draw_rectangle( first.xcoord(), first.ycoord(), x, y);
  w.set_buttons( save_but);
  return w;
}

template< class R >
leda_window&
read(leda_window& w, Iso_rectangle_2<R>& r)
{
  typedef typename R::RT    RT;
  double x,y;
  int key = 0;
#if ( __LEDA__ < 362 )
  w.state = 1;
#else
  w.set_state( 1);
#endif // __LEDA__ < ...
  leda_point first, second;
  leda_drawing_mode save = w.set_mode(leda_xor_mode);
  if ( !( w.read( first))) { w.set_mode( save); return w; }
  int save_but[8];
  w.std_buttons(save_but);
  key = w.read_mouse_rect( first.xcoord(), first.ycoord(), x, y);
  if ( key == MOUSE_BUTTON(3))
  {
      w.set_buttons( save_but);
      w.set_mode( save);
#if ( __LEDA__ < 362 )
      w.state = 0;
#else
      w.set_state( 0);
#endif // __LEDA__ < ...
      return w;
  }
  r = Iso_rectangle_2<R>(Point_2<R>( RT(first.xcoord()),
                                               RT(first.ycoord())),
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
leda_window&
operator<<(leda_window& w, const Bbox_2& b)
{
#if (__LEDA__ >= 400)
  leda_line_style style = w.set_line_style(leda_dotted);
#else
  line_style style = w.set_line_style(leda_dotted);
#endif
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

