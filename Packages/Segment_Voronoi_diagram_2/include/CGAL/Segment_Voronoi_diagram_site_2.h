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

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

CGAL_BEGIN_NAMESPACE

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class Gt>
class Segment_Voronoi_diagram_site_2 
{
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Point_2   Point_2;
  typedef typename Geom_traits::Segment_2 Segment_2;

protected:
  typedef typename Geom_traits::FT        FT;
  typedef typename Geom_traits::RT        RT;
  typedef Segment_Voronoi_diagram_site_2<Geom_traits>  Self;

public:
  Segment_Voronoi_diagram_site_2() : type_(0) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_site_2(const Point_2& p) {
    initialize_site(p);
  }

  // constructs segment site using the segment (p1,p2)
  Segment_Voronoi_diagram_site_2(const Point_2& p1, const Point_2& p2) {
    initialize_site(p1, p2);
  }

  // constructs point site using the point of intersection of the
  // segments (p1,p2) and (q1,q2)
  Segment_Voronoi_diagram_site_2(const Point_2& p1, const Point_2& p2,
				 const Point_2& q1, const Point_2& q2) {
    initialize_site(p1, p2, q1, q2);
  }

  // constructs segment site using the points of intersection of the
  // segment-pairs (p1,p2)-(q1,q2) and (p1,p2)-(q3,q4) as endpoints;
  // the segment (p1,p2) is a segment that supports the actual segment
  Segment_Voronoi_diagram_site_2(const Point_2& p1, const Point_2& p2,
				 const Point_2& q1, const Point_2& q2,
				 const Point_2& r1, const Point_2& r2) {
    initialize_site(p1, p2, q1, q2, r1, r2);
  }

  // constructs segment site using either the source or the target of
  // (p1,p2) (that depends on the boolean is_first_exact) and the
  // intersection of (p1,p2) with (q1,q2) as the other endpoint
  Segment_Voronoi_diagram_site_2(const Point_2& p1, const Point_2& p2,
				 const Point_2& q1, const Point_2& q2,
				 bool is_first_exact) {
    initialize_site(p1, p2, q1, q2, is_first_exact);
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

  const Point_2& point(unsigned int i) const
  {
    CGAL_precondition( i < 6 );
    if ( i == 0 ) { return p_[0]; }
    else if ( i == 1 ) {
      CGAL_precondition( is_segment() || !is_exact() );
      return p_[1];
    } else if ( i == 2 ) {
      CGAL_precondition( (is_point() && !is_exact()) ||
			 (is_segment() && !is_exact(0)) );
      return p_[2];
    } else if ( i == 3 ) {
      CGAL_precondition( (is_point() && !is_exact()) ||
			 (is_segment() && !is_exact(0)) );
      return p_[3];
    } else if ( i == 4 ) {
      CGAL_precondition( (is_point() && !is_exact()) ||
			 (is_segment() && !is_exact(1)) );
      return p_[4];
    } else {  // i == 5
      CGAL_precondition( (is_point() && !is_exact()) ||
			 (is_segment() && !is_exact(1)) );
      return p_[5];
    }
  }


  Point_2 point() const { 
    CGAL_precondition ( is_point() );
    if ( !is_exact() ) {
      return compute_intersection_point(p_[2], p_[3], p_[4], p_[5]);
    } else {
      return p_[0];
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

  Self supporting_site() const {
    CGAL_precondition( is_segment() );
    return Self(p_[0], p_[1]);
  }

  Self supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && i < 2);
    CGAL_precondition( !is_exact() );
    if ( i == 0 ) { return Self(p_[2], p_[3]); }
    return Self(p_[4], p_[5]);
  }

  Self crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_exact() );
    CGAL_precondition( i < 2 && !is_exact(i) );
    if ( i == 0 ) {
      return Self(p_[2], p_[3]);
    } else {
      return Self(p_[4], p_[5]);
    }
  }

  Self source_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(0) ) {
      return Self(p_[0]);
    } else {
      return Self(p_[0], p_[1], p_[2], p_[3]);
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(1) ) {
      return Self(p_[1]);
    } else {
      return Self(p_[0], p_[1], p_[4], p_[5]);
    }
  }

  Self opposite_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() ) {
      return Self(p_[1], p_[0]);
    }

    CGAL_assertion( !is_exact(0) || !is_exact(1) );

    if ( is_exact(0) && !is_exact(1) ) {
      return Self(p_[1], p_[0], p_[4], p_[5], false);
    } else if ( !is_exact(0) && is_exact(1) ) {
      return Self(p_[1], p_[0], p_[2], p_[3], true);
    } else {
      return Self(p_[1], p_[0], p_[4], p_[5], p_[2], p_[3]);
    }
  }

#if 1
  // MK::ERROR: at some point I want to remove these
  Segment_2 supporting_segment() const {
    CGAL_precondition( is_segment() );
    return Segment_2(p_[0], p_[1]);
  }

  Segment_2 supporting_segment(unsigned int i) const {
    CGAL_precondition( is_point() && !is_exact() && i < 2 );
    if ( i == 0 ) {
      return Segment_2(p_[2], p_[3]);
    } else {
      return Segment_2(p_[4], p_[5]);
    }
  }

  Segment_2 crossing_segment(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_exact() );
    CGAL_precondition( i < 2 && !is_exact(i) );
    if ( i == 0 ) {
      return Segment_2(p_[2], p_[3]);
    } else {
      return Segment_2(p_[4], p_[5]);
    }
  }
#endif

#ifdef USE_SET_METHODS
  void set_point(const Point_2& p) {
    initialize_site(p);
  }

  void set_segment(const Point_2& p1, const Point_2& p2) {
    initialize_site(p1, p2);
  }

  void set_point(const Point_2& p1, const Point_2& p2,
		 const Point_2& q1, const Point_2& q2) {
    initialize_site(p1, p2, q1, q2);
  }

  void set_segment(const Point_2& p1, const Point_2& p2,
		   const Point_2& q1, const Point_2& q2,
		   const Point_2& r1, const Point_2& r2) {
    initialize_site(p1, p2, q1, q2, r1, r2);
  }

  void set_segment(const Point_2& p1, const Point_2& p2,
		   const Point_2& q1, const Point_2& q2,
		   bool is_first_exact) {
    initialize_site(p1, p2, q1, q2, is_first_exact);
  }
#endif

protected:
  void initialize_site(const Point_2& p)
  {
    type_ = 1;
    p_[0] = p;
  }

  void initialize_site(const Point_2& p1, const Point_2& p2)
  {
    type_ = 2;
    p_[0] = p1;
    p_[1] = p2;
  }
  void initialize_site(const Point_2& p1, const Point_2& p2,
		       const Point_2& q1, const Point_2& q2)
  {
    // MK: Sort the segments s1 and s2 in lexicographical order so
    //     that the computation of the intersection point is always
    //     done in the same manner (?)
    type_ = 5;
    p_[2] = p1;
    p_[3] = p2;
    p_[4] = q1;
    p_[5] = q2;
  }


  void initialize_site(const Point_2& p1, const Point_2& p2,
		       const Point_2& q1, const Point_2& q2,
		       const Point_2& r1, const Point_2& r2)
  {
    type_ = 14;
    p_[0] = p1;
    p_[1] = p2;
    p_[2] = q1;
    p_[3] = q2;
    p_[4] = r1;
    p_[5] = r2;
  }

  void initialize_site(const Point_2& p1, const Point_2& p2,
		       const Point_2& q1, const Point_2& q2,
		       bool is_first_exact)
  {
    type_ = (is_first_exact ? 10 : 6);
    p_[0] = p1;
    p_[1] = p2;
    if ( is_first_exact ) {
      p_[4] = q1;
      p_[5] = q2;
    } else {
      p_[2] = q1;
      p_[3] = q2;
    }
  }

  Point_2 compute_source() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(0) ) {
      return p_[0];
    } else {
      return compute_intersection_point(p_[0], p_[1], p_[2], p_[3]);
    }
  }

  Point_2 compute_target() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(1) ) {
      return p_[1];
    } else {
      return compute_intersection_point(p_[0], p_[1], p_[4], p_[5]);
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
  Point_2 p_[6];
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

  char type;
  if (is >> type) {
    if (type == 'p') {
      Point_2 p;
      is >> p;
#ifdef USE_SET_METHODS
      t.set_point(p);
#else
      t = Site_2(p);
#endif
    } else if (type == 's') {
      Point_2 p1, p2;
      is >> p1 >> p2;
#ifdef USE_SET_METHODS
      t.set_segment(p1, p2);
#else
      t = Site_2(p1, p2);
#endif
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
