// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef SEGMENT_VORONOI_DIAGRAM_SITE_H
#define SEGMENT_VORONOI_DIAGRAM_SITE_H

#include <iostream>
#include <CGAL/assertions.h>

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
  Segment_Voronoi_diagram_site_2() : type_(0) {}

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
    if ( assign(p_, o) ) {
      initialize_site(p_);
      return;
    }

    Segment_2 s;
    if ( assign(s, o) ) {
      initialize_site(s);
      return;
    }

    type_ = 0;
  }

  bool is_defined() const { return type_; }
  bool is_point() const { return (type_ & 3) == 1; }
  bool is_segment() const { return (type_ & 3) == 2; }
  bool is_exact() const { return !(type_ & 12); }
  bool is_exact(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    if ( i == 0 ) { return !(type_ & 4); }
    return !(type_ & 8);
  }

  Point_2 point() const { 
    CGAL_precondition ( is_point() );
    if ( !is_exact() ) {
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
    if ( is_exact() || is_exact(0) ) {
      return Self(p_);
    } else {
      return Self(supporting_segment(), Segment_2(q1_, q2_));
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(1) ) {
      return Self(p2_);
    } else {
      return Self(supporting_segment(), Segment_2(q3_, q4_));
    }
  }

  Self opposite_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() ) {
      return Self( segment().opposite() );
    }

    Segment_2 supp = supporting_segment().opposite();

    CGAL_assertion( !is_exact(0) || !is_exact(1) );

    if ( is_exact(0) && !is_exact(1) ) {
      return Self(supp, crossing_segment(1), false);
    } else if ( !is_exact(0) && is_exact(1) ) {
      return Self(supp, crossing_segment(0), true);
    } else {
      return Self(supp, crossing_segment(1), crossing_segment(0));
    }
  }

  Segment_2 supporting_segment() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() ) {
      return segment();
    } else {
      return Segment_2(p_, p2_);
    }
  }

  Segment_2 supporting_segment(unsigned int i) const {
    CGAL_precondition( is_point() && !is_exact() && i < 2 );
    if ( i == 0 ) {
      return Segment_2(q1_, q2_);
    } else {
      return Segment_2(q3_, q4_);
    }
  }

  Segment_2 crossing_segment(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_exact() );
    CGAL_precondition( i < 2 && !is_exact(i) );
    if ( i == 0 ) {
      return Segment_2(q1_, q2_);
    } else {
      return Segment_2(q3_, q4_);
    }
  }

  void set_point(const Point_2& p) {
    initialize_site(p);
  }

  void set_segment(const Segment_2& s) {
    initialize_site(s);
  }

  void set_point(const Segment_2& s1, const Segment_2& s2) {
    initialize_site(s1, s2);
  }

  void set_segment(const Segment_2& support, const Segment_2& s1,
		   const Segment_2& s2) {
    initialize_site(support, s1, s2);
  }

  void set_segment(const Segment_2& support, const Segment_2& s,
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
    type_ = 1;
    p_ = p;
  }

  void initialize_site(const Segment_2& s)
  {
    type_ = 2;
    p_ = s.source();
    p2_ = s.target();
  }
  void initialize_site(const Segment_2& s1, const Segment_2& s2)
  {
    // MK: Sort the segments s1 and s2 in lexicographical order so
    //     that the computation of the intersection point is always
    //     done in the same manner (?)
    type_ = 5;
    q1_ = s1.source();
    q2_ = s1.target();
    q3_ = s2.source();
    q4_ = s2.target();
  }


  void initialize_site(const Segment_2& support, const Segment_2& s1,
		       const Segment_2& s2)
  {
    type_ = 14;
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
    type_ = (is_first_exact ? 10 : 6);
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
    if ( is_exact() || is_exact(0) ) {
      return p_;
    } else {
      return compute_intersection_point(p_, p2_, q1_, q2_);
    }
  }

  Point_2 compute_target() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(1) ) {
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
  char type_;
};

//-------------------------------------------------------------------------

template <class R>
std::ostream&
operator<<(std::ostream& os, 
	   const Segment_Voronoi_diagram_site_2<R>& s)
{
  if (!s.is_defined())
    return os << "u";
  if (s.is_point())
    return os << "p " << s.point ();
  return os << "s " << s.segment ();
}

template <class R>
std::istream &
operator>>(std::istream &is,
	   Segment_Voronoi_diagram_site_2<R>& t)
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
