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

protected:
  typedef typename R::FT        FT;
  typedef typename R::RT        RT;
  typedef Segment_Voronoi_diagram_site_2<Rep>  Self;

public:
  Segment_Voronoi_diagram_site_2() : defined_(false) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_site_2(const Point_2 &p) {
    initialize_site(p);
  }

  // constructs segment site using input segment
  Segment_Voronoi_diagram_site_2(const Segment_2 &s) {
    initialize_site(s);
  }

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

  bool is_defined () const { return defined_; }
  bool is_point () const { return defined_ && point_; }
  bool is_segment () const { return defined_ && !point_; }
  bool is_exact() const { return input_; }
  bool is_exact(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    return is_exact_[i];
  }

  Point_2 point() const { 
    CGAL_precondition ( is_point() );
    if ( !input_ ) {
      return compute_intersection_point(q1_, q2_, q3_, q4_);
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

  Self source_site() const {
    CGAL_precondition( is_segment() );
    if ( input_ || is_exact_[0] ) {
      return Self(p_);
    } else {
      return Self(supporting_segment(), Segment_2(q1_, q2_));
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( input_ || is_exact_[1] ) {
      return Self(p2_);
    } else {
      return Self(supporting_segment(), Segment_2(q3_, q4_));
    }
  }

  Segment_2 supporting_segment() const {
    CGAL_precondition( is_segment() );
    if ( input_ ) {
      return segment();
    } else {
      return Segment_2(p_, p2_);
    }
  }

  Segment_2 supporting_segment(unsigned int i) const {
    CGAL_precondition( is_point() && !input_ && i < 2 );
    if ( i == 0 ) {
      return Segment_2(q1_, q2_);
    } else {
      return Segment_2(q3_, q4_);
    }
  }

  Segment_2 crossing_segment(unsigned int i) const {
    CGAL_precondition( is_segment() && !input_ );
    CGAL_precondition( i < 2 && !is_exact_[i] );
    if ( i == 0 ) {
      return Segment_2(q1_, q2_);
    } else {
      return Segment_2(q3_, q4_);
    }
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
  void initialize_site(const Point_2& p)
  {
    defined_ = true;
    point_ = true;
    input_ = true;
    is_exact_[0] = is_exact_[1] = false;
    p_ = p;
  }

  void initialize_site(const Segment_2& s)
  {
    defined_ = true;
    point_ = false;
    input_ = true;
    is_exact_[0] = is_exact_[1] = true;
    p_ = s.source();
    p2_ = s.target();
  }
  void initialize_site(const Segment_2& s1, const Segment_2& s2)
  {
    // MK: Sort the segments s1 and s2 in lexicographical order so
    //     that the computation of the intersection point is always
    //     done in the same manner (?)
    defined_ = true;
    point_ = true;
    input_ = false;
    is_exact_[0] = is_exact_[1] = false;
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
    is_exact_[0] = is_exact_[1] = false;
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
    is_exact_[0] = is_first_exact;
    is_exact_[1] = !is_first_exact;
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
    CGAL_precondition( is_segment() );
    if ( input_ || is_exact_[0] ) {
      return p_;
    } else {
      return compute_intersection_point(p_, p2_, q1_, q2_);
    }
  }

  Point_2 compute_target() const {
    CGAL_precondition( is_segment() );
    if ( input_ || is_exact_[1] ) {
      return p2_;
    } else {
      return compute_intersection_point(p_, p2_, q3_, q4_);
    }
  }

  // computes the point of intersection of the segments p1p2 and p3p4
  static Point_2
  compute_intersection_point(const Point_2& p1, const Point_2& p2,
			     const Point_2& p3, const Point_2& p4)
  {
    RT x1 = p1.x(), y1 = p1.y();
    RT x2 = p2.x(), y2 = p2.y();
    RT x3 = p3.x(), y3 = p3.y();
    RT x4 = p4.x(), y4 = p4.y();

    RT D = det2x2_by_formula(x2 - x1, x4 - x3, y2 - y1, y4 - y3);
    RT Dt = det2x2_by_formula(x3 - x1, x4 - x3, y3 - y1, y4 - y3);

    RT t = Dt / D;

    return Point_2(x1 + (x2 - x1) * t, y1 + (y2 - y1) * t);
  }

protected:
  Point_2 p_;
  Point_2 p2_;
  Point_2 q1_, q2_, q3_, q4_;
  bool defined_;
  bool point_;
  bool input_;
  bool is_exact_[2];
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
