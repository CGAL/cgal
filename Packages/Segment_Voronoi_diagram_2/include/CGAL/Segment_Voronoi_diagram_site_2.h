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

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class R_>
class Segment_Voronoi_diagram_site_2 
{
public:
  typedef R_ R;
  typedef R  Rep;
  typedef typename R::Point_2   Point_2;
  typedef typename R::Segment_2 Segment_2;

  Segment_Voronoi_diagram_site_2() : defined_(false) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_site_2(const Point_2 &p)
    : p_(p), defined_(true), point_(true), input_(true)  {}

  // constructs segment site using input segment
  Segment_Voronoi_diagram_site_2(const Segment_2 &s)
    : p_(s.source()), p2_(s.target()), defined_(true), point_(false),
      input_(true), is_exact1_(true), is_exact2_(true) {}

  // constructs point site using point of intersection
  Segment_Voronoi_diagram_site_2(const Segment_2& s1,
				 const Segment_2& s2) {
    initialize_site(s1, s2);
  }

  // constructs segment site using points of intersection of support
  // with s1 and support with s2 as endpoints
  Segment_Voronoi_diagram_site_2(const Segment_2& support,
				 const Segment_2& s1,
				 const Segment_2& s2) {
    initialize_site(support, s1, s2);
  }

  // constructs segment site using either the source or the target of
  // support (that depends on the boolean is_first_exact) and the
  // intersection of support with s as the other endpoint
  Segment_Voronoi_diagram_site_2(const Segment_2& support,
				 const Segment_2& s,
				 bool is_first_exact) {
    initialize_site(support, s, is_first_exact);
  }

  Segment_Voronoi_diagram_site_2(const Object &o) {
    if ( assign(_p, o) ) {
      defined_ = true;
      point_ = true;
      input_ = true;
      return;
    }

    Segment_2 s;
    if ( assign(s, o) ) {
      p_ = s.source();
      p2_ = s.target();
      defined_ = true;
      point_ = false;
      input_ = true;
      return;
    }

    defined_ = false;
  }

  //***********************************************************************

  bool is_defined () const { return defined_; }
  bool is_point () const { return defined_ && point_; }
  bool is_segment () const { return defined_ && !point_; }
  bool is_exact() const { return input_; }

  Point_2 point() const { 
    CGAL_precondition ( is_point() );
    if ( !input_ ) {
      return compute_intersection_point();
    } else {
      return p_;
    }
  }
  Segment_2 segment() const {
    CGAL_precondition ( is_segment() ); 
    return Segment_2( source(), target() );
  }

  Point_2 source() const {
    CGAL_precondition ( is_segment() ); 
    return compute_source();
  }

  Point_2 target() const {
    CGAL_precondition ( is_segment() ); 
    return compute_target();
  }

  void set_point(const Point_2& p) {
    p_ = p;
    defined_ = true;
    point_ = true;
    input_ = true;
  }

  void set_segment(const Segment_2& s) {
    p_ = s.source();
    p2_ = s.target();
    defined_ = true;
    point_ = false;
    input_ = true;
  }

  void set_point(const Segment_2& s1, Segment_2& s2) {
    initialize_site(s1, s2);
  }

  void set_segment(const Segment_2& support, const Segment_2& s1,
		   Segment_2& s2) {
    initialize_site(support, s1, s2);
  }

  void set_point(const Segment_2& support, Segment_2& s,
		 bool is_first_exact) {
    initialize_site(support, s, is_first_exact);
  }

  std::ostream& write(std::ostream& os)
  {
    return os << (*this);
  }

protected:
  void initialize_site(const Segment_2& s1, const Segment_2& s2)
  {
    defined_ = true;
    point_ = true;
    input_ = false;
    q1_ = s1.source();
    q2_ = s1.target();
    q3_ = s2.source();
    q4_ = s2.target();
  }


  void initialize_site(const Segment_2& support, const Segment_2& s1,
		       const Segment_2& s2)
  {
    defined_ = true;
    point_ = false;
    input_ = false;
    is_exact1_ = false;
    is_exact2_ = false;
    p_ = support.source();
    p2_ = support.target();
    q1_ = s1.source();
    q2_ = s1.target();
    q3_ = s2.source();
    q4_ = s2.target();
  }

  void initialize_site(const Segment_2& support, const Segment_2& s,
		       bool is_first_exact)
  {
    defined_ = true;
    point_ = false;
    input_ = false;
    is_exact1_ = is_first_exact;
    is_exact2_ = !is_first_exact;
    p_ = support.source();
    p2_ = support.target();
    if ( is_first_exact ) {
      q3_ = s.source();
      q4_ = s.target();
    } else {
      q1_ = s.source();
      q2_ = s.target();
    }
  }

  Point_2 compute_source() const {
    if ( input_ || is_exact1_ ) {
      return p_;
    } else {
      return compute_intersection_point(0);
    }
  }

  Point_2 compute_target() const {
    if ( input_ || is_exact2_ ) {
      return p2_;
    } else {
      return compute_intersection_point(1);
    }
  }

  // MK: the following two methods have to be filled in
  Point_2 compute_intersection_point() const
  {
    CGAL_precondition( point_ && !input_ );
    return Point_2();
  }

  Point_2 compute_intersection_point(int i) const
  {
    CGAL_precondition( !point_ && !input_ );
    CGAL_precondition( (i == 0 && !is_exact1_) ||
		       (i == 1 && !is_exact2_) );
    return Point_2();
  }

protected:
  Point_2 p_;
  Point_2 p2_;
  Point_2 q1_, q2_, q3_, q4_;
  bool defined_;
  bool point_;
  bool input_;
  bool is_exact1_, is_exact2_;
};

//-------------------------------------------------------------------------

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


CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_SITE_H
