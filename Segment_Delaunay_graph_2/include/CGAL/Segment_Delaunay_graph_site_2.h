// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_SITE_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_SITE_H

#include <iostream>
#include <CGAL/assertions.h>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_2/Constructions_C2.h>

namespace CGAL {

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class Gt>
class Segment_Delaunay_graph_site_2 
{
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Point_2   Point_2;
  typedef typename Geom_traits::Segment_2 Segment_2;

protected:
  typedef typename Geom_traits::FT        FT;
  typedef typename Geom_traits::RT        RT;
  typedef Segment_Delaunay_graph_site_2<Geom_traits>  Self;

public:
  Segment_Delaunay_graph_site_2() : type_(0) {}

  static Self construct_site_2(const Point_2& p) {
    Self t;
    t.initialize_site(p);
    return t;
  }

  // constructs segment site using the segment (p1,p2)
  static Self construct_site_2(const Point_2& p1, const Point_2& p2) {
    Self t;
    t.initialize_site(p1, p2);
    return t;
  }

  // constructs point site using the point of intersection of the
  // segments (p1,p2) and (q1,q2)
  static Self construct_site_2(const Point_2& p1, const Point_2& p2,
			       const Point_2& q1, const Point_2& q2) {
    Self t;
    t.initialize_site(p1, p2, q1, q2);
    return t;
  }

  // constructs segment site using the points of intersection of the
  // segment-pairs (p1,p2)-(q1,q2) and (p1,p2)-(q3,q4) as endpoints;
  // the segment (p1,p2) is a segment that supports the actual segment
  static Self construct_site_2(const Point_2& p1, const Point_2& p2,
			       const Point_2& q1, const Point_2& q2,
			       const Point_2& r1, const Point_2& r2) {
    Self t;
    t.initialize_site(p1, p2, q1, q2, r1, r2);
    return t;
  }

  // constructs segment site using either the source or the target of
  // (p1,p2) (that depends on the boolean is_first_exact) and the
  // intersection of (p1,p2) with (q1,q2) as the other endpoint
  static Self construct_site_2(const Point_2& p1, const Point_2& p2,
			       const Point_2& q1, const Point_2& q2,
			       bool is_first_exact) {
    Self t;
    t.initialize_site(p1, p2, q1, q2, is_first_exact); 
    return t;
 }

public:
  bool is_defined() const { return type_ != 0; }
  bool is_point() const { return (type_ & 3) == 1; }
  bool is_segment() const { return (type_ & 3) == 2; }
  bool is_input() const { return !(type_ & 12); }
  bool is_input(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    if ( i == 0 ) { return !(type_ & 4); }
    return !(type_ & 8);
  }

  const Point_2& source_of_supporting_site() const {
    CGAL_precondition( is_segment() );
    return p_[0];
  }

  const Point_2& target_of_supporting_site() const {
    CGAL_precondition( is_segment() );
    return p_[1];
  }

  const Point_2& source_of_supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() );
    return (i == 0) ? p_[2] : p_[4];
  }

  const Point_2& target_of_supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() );
    return (i == 0) ? p_[3] : p_[5];
  }

  const Point_2& source_of_crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input(i) );
    return (i == 0) ? p_[2] : p_[4];
  }

  const Point_2& target_of_crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input(i) );
    return (i == 0) ? p_[3] : p_[5];
  }

  Point_2 point() const { 
    CGAL_precondition ( is_point() );
    if ( !is_input() ) {
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
    return construct_site_2(p_[0], p_[1]);
  }

  Self supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && i < 2);
    CGAL_precondition( !is_input() );
    if ( i == 0 ) { return construct_site_2(p_[2], p_[3]); }
    return construct_site_2(p_[4], p_[5]);
  }

  Self crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input() );
    CGAL_precondition( i < 2 && !is_input(i) );
    if ( i == 0 ) {
      return construct_site_2(p_[2], p_[3]);
    } else {
      return construct_site_2(p_[4], p_[5]);
    }
  }

  Self source_site() const {
    CGAL_precondition( is_segment() );
    if ( is_input() || is_input(0) ) {
      return construct_site_2(p_[0]);
    } else {
      return construct_site_2(p_[0], p_[1], p_[2], p_[3]);
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( is_input() || is_input(1) ) {
      return construct_site_2(p_[1]);
    } else {
      return construct_site_2(p_[0], p_[1], p_[4], p_[5]);
    }
  }

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
    if ( is_input() || is_input(0) ) {
      return p_[0];
    } else {
      return compute_intersection_point(p_[0], p_[1], p_[2], p_[3]);
    }
  }

  Point_2 compute_target() const {
    CGAL_precondition( is_segment() );
    if ( is_input() || is_input(1) ) {
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

    RT D = determinant(x2 - x1, x4 - x3, y2 - y1, y4 - y3);
    RT Dt = determinant(x3 - x1, x4 - x3, y3 - y1, y4 - y3);

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
	   const Segment_Delaunay_graph_site_2<R>& s)
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
	   Segment_Delaunay_graph_site_2<R>& t)
{
  typedef Segment_Delaunay_graph_site_2<R>   Site_2;
  typedef typename Site_2::Point_2           Point_2;

  char type;
  if (is >> type) {
    if (type == 'p') {
      Point_2 p;
      is >> p;
      t = Site_2::construct_site_2(p);
    } else if (type == 's') {
      Point_2 p1, p2;
      is >> p1 >> p2;
      t = Site_2::construct_site_2(p1, p2);
    }
  }
  return is;
}

template < class R, class Stream >
Stream&
operator<<(Stream& str, Segment_Delaunay_graph_site_2<R>& t)
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


} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_SITE_H
