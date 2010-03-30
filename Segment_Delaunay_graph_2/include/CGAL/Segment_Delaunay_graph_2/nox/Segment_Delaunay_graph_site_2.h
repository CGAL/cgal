// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France) and
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_NOX_SITE_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_NOX_SITE_H

#include <iostream>
#include <CGAL/assertions.h>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>


CGAL_BEGIN_NAMESPACE

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
  typedef typename Geom_traits::Line_2    Line_2;
  typedef Segment_Delaunay_graph_site_2<Geom_traits>  Self;

public:
  Segment_Delaunay_graph_site_2() : type_(0) {}

#if 0
  Segment_Delaunay_graph_site_2(const Point_2& p)
    : type_(1), p1_(p) {}
#ifdef CGAL_SDG_SORT_POINTS_IN_SITE2
  Segment_Delaunay_graph_site_2(const Point_2& p1, const Point_2& p2)
    : type_(2)
  {
    Comparison_result rs = CGAL::compare(p1.x(), p2.x());
    if ( rs == SMALLER ) {
      p1_ = p1;
      p2_ = p2;
    } else if ( rs == LARGER ) {
      p1_ = p2;
      p2_ = p1;
    } else {  // rs == EQUAL
      rs = CGAL::compare(p1.y(), p2.y());
      CGAL_assertion( rs != EQUAL );
      if ( rs == SMALLER ) {
	p1_ = p1;
	p2_ = p2;
      } else {
	p1_ = p2;
	p2_ = p1;
      }
    }
  }
#else
  Segment_Delaunay_graph_site_2(const Point_2& p1, const Point_2& p2)
    : type_(2), p1_(p1), p2_(p2) {}
#endif
#endif

public:
  inline bool is_defined() const { return type_ != 0; }
  inline bool is_point()   const { return type_ == 1; }
  inline bool is_segment() const { return type_ == 2; }

  inline const Point_2& point() const { 
    CGAL_precondition ( is_point() );
    return p1_;
  }

  inline Segment_2 segment() const {
    CGAL_precondition ( is_segment() ); 
    return Segment_2( p1_, p2_ );
  }

  inline const Point_2& source() const {
    CGAL_precondition ( is_segment() ); 
    return p1_;
  }

  inline const Point_2& target() const {
    CGAL_precondition ( is_segment() ); 
    return p2_;
  }

#if 1
  static Self construct_site_2(const Point_2& p) {
    Self t;
    t.initialize_site(p);
    return t;
  }

  static Self construct_site_2(const Point_2& p1, const Point_2& p2) {
    Self t;
    t.initialize_site(p1, p2);
    return t;
  }
#else
  static Self construct_site_2(const Point_2& p) { return Self(p); }
  static Self construct_site_2(const Point_2& p1, const Point_2& p2) {
    return Self(p1, p2);
  }
#endif

protected:
  void initialize_site(const Point_2& p)
  {
    type_ = 1;
    p1_ = p;
  }

  void initialize_site(const Point_2& p1, const Point_2& p2)
  {
    type_ = 2;
#ifdef CGAL_SDG_SORT_POINTS_IN_SITE2
    Comparison_result rs = CGAL::compare(p1.x(), p2.x());
    if ( rs == SMALLER ) {
      p1_ = p1;
      p2_ = p2;
    } else if ( rs == LARGER ) {
      p1_ = p2;
      p2_ = p1;
    } else {  // rs == EQUAL
      rs = CGAL::compare(p1.y(), p2.y());
      CGAL_assertion( rs != EQUAL );
      if ( rs == SMALLER ) {
	p1_ = p1;
	p2_ = p2;
      } else {
	p1_ = p2;
	p2_ = p1;
      }
    }
#else
    p1_ = p1;
    p2_ = p2;
#endif
  }


protected:
  char type_;
  Point_2 p1_, p2_;
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
    return os << "p " << s.point();
  return os << "s " << s.segment();
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
      t = Site_2(p);
    } else if (type == 's') {
      Point_2 p1, p2;
      is >> p1 >> p2;
      t = Site_2(p1, p2);
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
      str << "s " << t.source() << "  "	<< t.target();
    }
  }

  return str;
}


CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_NOX_SITE_H
