// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Segment_Voronoi_diagram_site_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef SEGMENT_VORONOI_DIAGRAM_SITE_H
#define SEGMENT_VORONOI_DIAGRAM_SITE_H

#include <iostream>
#include <CGAL/assertions.h>

#ifdef VORONOI_WANT_WINDOW
#  include <CGAL/IO/Window_stream.h>
#endif

CGAL_BEGIN_NAMESPACE

  /** A Site is either a point or a segment (if defined)
   */

template <class R_>
class Segment_Voronoi_diagram_site_2 
{
public:
  typedef R_ R;
  typedef R  Rep;
  typedef typename R::Point_2   Point_2;
  typedef typename R::Segment_2 Segment_2;
  typedef Segment_Voronoi_diagram_site_2<R>  Site_2;

  Segment_Voronoi_diagram_site_2 () : defined_ (false) {}
  ~Segment_Voronoi_diagram_site_2 () { }

  Segment_Voronoi_diagram_site_2 (const Point_2 &p)
    : _p(p), defined_ (true), point_ (true)  {}
  
  Segment_Voronoi_diagram_site_2 (const Segment_2 &s)
    : _p(s.source()), _p2(s.target()), defined_ (true), point_ (false) { }

  Segment_Voronoi_diagram_site_2 (const Object &o) {
    if ( assign(_p, o) ) {
      defined_ = true;
      point_ = true;
      return;
    }

    Segment_2 s;
    if ( assign(s, o) ) {
      //  *this = s;
      _p = s.source();
      _p2 = s.target();
      defined_ = true;
      point_ = false;
      return;
    }

    defined_ = false;
  }

  inline bool is_defined () const { return defined_; }
  inline bool is_point () const { return defined_ && point_; }
  inline bool is_segment () const { return defined_ && !point_; }

  inline const Point_2& point() const { 
    CGAL_precondition ( is_point() );
    return _p;
  }
  inline const Segment_2 segment() const {
    CGAL_precondition ( is_segment() ); 
    return Segment_2(_p,_p2);
  }

  inline const Point_2& source() const {
    CGAL_precondition ( is_segment() ); 
    return _p;
  }

  inline const Point_2& target() const {
    CGAL_precondition ( is_segment() ); 
    return _p2;
  }

#if 0
  inline Point_2& point() { 
    CGAL_precondition ( is_point() );
    return _p;
  }
  inline Segment_2& segment() {
    CGAL_precondition ( is_segment() ); 
    return _s;
  }
#endif

  void set_point(const Point_2& p) {
    _p = p;
    defined_ = true;
    point_ = true;
  }

  void set_segment(const Segment_2& s) {
    _p = s.source();
    _p2 = s.target();
    defined_ = true;
    point_ = false;
  }

  inline std::ostream& write(std::ostream& os)
  {
    return os << (*this);
  }

protected:
  Point_2 _p;
  Point_2 _p2;
  bool defined_;
  bool point_;
};


template <class R>
std::ostream&
operator<< (std::ostream& os, const Segment_Voronoi_diagram_site_2<R>& s)
{
  if (!s.is_defined())
    return os << "u";
  if (s.is_point())
    return os << "p " << s.point ();
  return os << "s " << s.segment ();
}

template <class R>
std::istream &
operator>> (std::istream &is, Segment_Voronoi_diagram_site_2<R>& t)
{
  typedef Segment_Voronoi_diagram_site_2<R>   Site_2;
  typedef typename Site_2::Point_2            Point_2;
  typedef typename Site_2::Segment_2          Segment_2;

  char type;
  if (is >> type) {
    if (type == 'p') {
      Point_2 p;
      is >> p;
      t.set_point(p);
    } else if (type == 's') {
      Segment_2 s;
      is >> s;
      t.set_segment(s);
    }
  }
  return is;
}

template < class R, class Stream >
Stream&
operator<<(Stream& str, Segment_Voronoi_diagram_site_2<R>& t)
{
  if ( t.is_defined() ) {
    if ( t.is_point() ) {
      str << "p " << t.point();
    } else {
      str << "s " << t.segment().source() << "  "
	  << t.segment().target();
    }
  }

  return str;
}

#ifdef VORONOI_WANT_WINDOW

#ifdef CHR_BIG_POINT
#  define LittleCircle(c,rad) (Circle_2<R> (c, (rad)*(rad)))
#  include <CGAL/Circle_2.h>
#endif

template <class R>
Window_stream &
operator<< (Window_stream &ws, const Site_2<R> &s)
{
  if (!s.is_defined())
    return ws;
  if (s.is_point()) {
#ifdef CHR_BIG_POINT
    return ws << LittleCircle (s.point (), 0.01)
	      << LittleCircle (s.point (), 0.02);
#else
    return ws << s.point ();
#endif
  }
  return ws << s.segment ();
}


template <class R>
Window_stream &
operator>> (Window_stream &ws, Site_2<R> &s)
{
  Segment_2<R> seg;

  if (! (ws >> seg)) return ws;
  
  if (seg.is_degenerate ()) {
    s = Site_2<R> (seg.start ());
    ws << s;
  } else
    s = Site_2<R> (seg);

  return ws;
}

#endif // VORONOI_WANT_WINDOW


CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_SITE_H
