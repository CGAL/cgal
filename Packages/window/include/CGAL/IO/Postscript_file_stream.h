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
// release       : 
// release_date  : 
// 
// file          : CGAL/IO/Postscript_file_stream.h
// package       : window
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POSTSCRIPT_FILE_STREAM_H
#define CGAL_POSTSCRIPT_FILE_STREAM_H

#include <CGAL/IO/Color.h>
#ifdef LEDA_PS_FILE_H
# error Internal CGAL error: <LEDA/ps_file.h> should not have been included yet
#else
# include <LEDA/basic.h>
# include <LEDA/list.h>
# include <LEDA/string.h>
# include <LEDA/stream.h>
# include <LEDA/point.h>
# include <LEDA/segment.h>
# include <LEDA/line.h>
# include <LEDA/circle.h>
# include <LEDA/polygon.h>
# if (__LEDA__ >= 400)
#  include <LEDA/rectangle.h>
# endif // 400
# include <LEDA/color.h>
# include <LEDA/window.h>
# if (__LEDA__ >= 400)
#  include <LEDA/rat_window.h>
#  include <LEDA/geo_graph.h>
# endif // 400
# define private protected
# include <LEDA/ps_file.h>
# undef private
#endif // LEDA_PS_FILE_H
#include <CGAL/IO/esprit_logo.xpm>

CGAL_BEGIN_NAMESPACE


typedef ::ps_file   leda_ps_file;

class Postscript_file_stream : public leda_ps_file
{
 public:
   Postscript_file_stream(double w,double h, leda_string name="CGAL_unnamed.ps")
   : leda_ps_file(w/40.0, h/40.0, name)
   {}

   Postscript_file_stream(leda_string name="CGAL_unnamed.ps")
    : leda_ps_file( name )
   { set_draw_bb(false); }

   void display() {}
   void display(int, int) {}
   int read_mouse(double& , double& ) {return 1;}
   leda_color set_fg_color(leda_color c) { return set_color(c); }
   void set_font(const leda_string& ls) { set_text_font(ls); }
   void change_rgb(const Color&);
   bool is_in_colormode();
};


inline
void
Postscript_file_stream::change_rgb(const Color& c)
{
  file << double(c.r())/255.0 << " "
       << double(c.g())/255.0 << " "
       << double(c.b())/255.0 << " srgb\n";
}


inline
bool
Postscript_file_stream::is_in_colormode()
{ return (outputmode==colored_mode); }


inline
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Color& c)
{
  assert( w.is_in_colormode() );
  w.change_rgb(c);
  return w;
}


inline
void
cgalize(Postscript_file_stream& w)
{
  w.set_line_width( 1);
}


inline
Postscript_file_stream*
create_demo_postscript_file_stream( float w = 512.0, float h = 512.0,
                                    const char* str = "CGAL_unnamed.ps",
                                    double x_extension = 1.0)
{
  Postscript_file_stream* Wptr = new Postscript_file_stream( w, h, str);
  cgalize( *Wptr);
  double y_extension = x_extension * h / w;
  Wptr->init(-x_extension, x_extension, -y_extension);
  return Wptr;
}


inline
Postscript_file_stream*
create_and_display_demo_postscript_file_stream(float w = 512.0, 
                    float h = 512.0, const char* str = "CGAL_unnamed.ps",
                    double x_extension = 1.0)
{
  Postscript_file_stream* Wptr = new Postscript_file_stream( w, h, str);
  cgalize( *Wptr);
  double y_extension = x_extension * h / w;
  Wptr->init(-x_extension, x_extension, -y_extension);
  return Wptr;
}


#endif // CGAL_POSTSCRIPT_FILE_STREAM_H

//  Each of the following operators is individually
//  protected against multiple inclusion.

#ifdef CGAL_POINT_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_POINT_2
#define CGAL_POSTSCRIPT_FILE_STREAM_POINT_2
template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Point_2<R>& p)
{
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  w.draw_point(x,y);

  return w;
}

template< class R >
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Point_2<R>& p)
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
Postscript_file_stream&
read(Postscript_file_stream& w, Point_2<R>& p)
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
read_mouse_plus(Postscript_file_stream& w, Point_2<R>& p, int& button)
{
  typedef typename R::RT RT;
  double x, y;
  button = w.read_mouse(x,y);
  w.draw_point(x,y);

  p = Point_2<R>(RT(x), RT(y));
}

#endif // CGAL_POSTSCRIPT_FILE_STREAM_POINT_2
#endif // CGAL_POINT_2_H


#ifdef CGAL_SEGMENT_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_SEGMENT_2
#define CGAL_POSTSCRIPT_FILE_STREAM_SEGMENT_2
template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Segment_2<R>& s)
{
  w.draw_segment(CGAL::to_double(s.source().x()),
                 CGAL::to_double(s.source().y()),
                 CGAL::to_double(s.target().x()),
                 CGAL::to_double(s.target().y()));
  return w;
}

template< class R >
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Segment_2<R>& s)
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
Postscript_file_stream&
read(Postscript_file_stream& w, Segment_2<R>& s)
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_SEGMENT_2
#endif // CGAL_SEGMENT_2_H


#ifdef CGAL_LINE_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_LINE_2
#define CGAL_POSTSCRIPT_FILE_STREAM_LINE_2

template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Line_2<R>& l)
{
  Point_2<R> p1 = l.point(),
                  p2 = p1 + l.direction().vector();
  w.draw_line(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),
              CGAL::to_double(p2.x()), CGAL::to_double(p2.y()));
  return w;
}

template< class R >
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Line_2<R>& l)
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
Postscript_file_stream&
read(Postscript_file_stream& w, Line_2<R>& l)
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_LINE_2
#endif //CGAL_LINE_2_H

#ifdef CGAL_RAY_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_RAY_2
#define CGAL_POSTSCRIPT_FILE_STREAM_RAY_2
template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Ray_2<R>& r)
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
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Ray_2<R>& r)
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
Postscript_file_stream&
read(Postscript_file_stream& w, Ray_2<R>& r)
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_RAY_2
#endif //CGAL_RAY_2_H

#ifdef CGAL_TRIANGLE_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_TRIANGLE_2
#define CGAL_POSTSCRIPT_FILE_STREAM_TRIANGLE_2
template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Triangle_2<R>& t)
{
  double x0 = CGAL::to_double(t.vertex(0).x()),
         y0 = CGAL::to_double(t.vertex(0).y()),
         x1 = CGAL::to_double(t.vertex(1).x()),
         y1 = CGAL::to_double(t.vertex(1).y()),
         x2 = CGAL::to_double(t.vertex(2).x()),
         y2 = CGAL::to_double(t.vertex(2).y());
  w.draw_segment(x0, y0, x1, y1);
  w.draw_segment(x1, y1, x2, y2);
  w.draw_segment(x2, y2, x0, y0);
  return w;
}

template< class R >
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Triangle_2<R>& t)
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
  w.draw_segment(x0,y0, x1, y1);
  w.draw_segment(x1,y1, x2, y2);
  w.draw_segment(x2,y2, x0, y0);

  t = Triangle_2<R>(Point_2<R>( RT(x0), RT(y0)),
                         Point_2<R>( RT(x1), RT(y1)),
                         Point_2<R>( RT(x2), RT(y2)));
  return w;
}

template< class R >
Postscript_file_stream&
read(Postscript_file_stream& w, Triangle_2<R>& t)
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_TRIANGLE_2
#endif // CGAL_TRIANGLE_2_H

#ifdef CGAL_CIRCLE_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_CIRCLE_2
#define CGAL_POSTSCRIPT_FILE_STREAM_CIRCLE_2
template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Circle_2<R>& c)
{
  double cx = CGAL::to_double(c.center().x()),
         cy = CGAL::to_double(c.center().y()),
         r = CGAL::to_double(c.squared_radius());
  w.draw_circle(cx, cy , sqrt(r));
  return w;
}

template< class R >
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Circle_2<R>& c)
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
  w.draw_circle(cx, cy , sqrt(sqr));
  c = Circle_2<R>(center, RT(sqr));
  return w;
}

template< class R >
Postscript_file_stream&
read(Postscript_file_stream& w, Circle_2<R>& c)
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_CIRCLE_2
#endif // CGAL_CIRCLE_2_H

#ifdef CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_ISO_RECTANGLE_2
#define CGAL_POSTSCRIPT_FILE_STREAM_ISO_RECTANGLE_2
template< class R >
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Iso_rectangle_2<R>& r)
{
  double xmin = CGAL::to_double(r.min().x()),
         ymin = CGAL::to_double(r.min().y()),
         xmax = CGAL::to_double(r.max().x()),
         ymax = CGAL::to_double(r.max().y());
  w.draw_segment(xmin, ymin, xmax, ymin);
  w.draw_segment(xmax, ymin, xmax, ymax);
  w.draw_segment(xmax, ymax, xmin, ymax);
  w.draw_segment(xmin, ymax, xmin, ymin);
  return w;
}

template< class R >
Postscript_file_stream&
operator>>(Postscript_file_stream& w, Iso_rectangle_2<R>& r)
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
  w.draw_rectangle( first.xcoord(), first.ycoord(), x, y);
  w.set_buttons( save_but);
  return w;
}

template< class R >
Postscript_file_stream&
read(Postscript_file_stream& w, Iso_rectangle_2<R>& r)
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_ISO_RECTANGLE_2
#endif // CGAL_ISO_RECTANGLE_2_H

#ifdef CGAL_BBOX_2_H
#ifndef CGAL_POSTSCRIPT_FILE_STREAM_BBOX_2
#define CGAL_POSTSCRIPT_FILE_STREAM_BBOX_2
inline
Postscript_file_stream&
operator<<(Postscript_file_stream& w, const Bbox_2& b)
{
#if (__LEDA__ >= 400)
  leda_line_style style = w.set_line_style(leda_dotted);
#else
  line_style style = w.set_line_style(leda_dotted);
#endif // 400
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
#endif // CGAL_POSTSCRIPT_FILE_STREAM_BBOX_2
#endif // CGAL_BBOX_2_H

CGAL_END_NAMESPACE

