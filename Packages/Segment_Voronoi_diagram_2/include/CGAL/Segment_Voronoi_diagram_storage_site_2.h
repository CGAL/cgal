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



#ifndef SEGMENT_VORONOI_DIAGRAM_STORAGE_SITE_H
#define SEGMENT_VORONOI_DIAGRAM_STORAGE_SITE_H

#include <iostream>
#include <CGAL/assertions.h>


CGAL_BEGIN_NAMESPACE

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class R_, class H_>
class Segment_Voronoi_diagram_storage_site_2 
{
public:
  typedef R_ R;
  typedef R  Rep;
  typedef H_ Handle;
  typedef typename R::Point_2   Point_2;
  typedef typename R::Segment_2 Segment_2;
  typedef typename R::Site_2    Site_2;

  typedef std::pair<Handle,Handle>  Handle_pair;
protected:
  typedef typename R::FT        FT;
  typedef typename R::RT        RT;
  typedef Segment_Voronoi_diagram_storage_site_2<Rep,Handle>  Self;

public:
  Segment_Voronoi_diagram_storage_site_2() : type_(0) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_storage_site_2(const Handle &h) {
    initialize_site(h);
  }

  // constructs segment site using input segment
  Segment_Voronoi_diagram_storage_site_2(const Handle_pair &hp) {
    initialize_site(hp);
  }

  // constructs point site using point of intersection
  Segment_Voronoi_diagram_storage_site_2(const Handle_pair& hp1,
					 const Handle_pair& hp2) {
    initialize_site(hp1, hp2);
  }

  // constructs segment site using points of intersection of support
  // with s1 and support with s2 as endpoints
  Segment_Voronoi_diagram_storage_site_2(const Handle_pair& hsupport,
					 const Handle_pair& hp1,
					 const Handle_pair& hp2) {
    initialize_site(hsupport, hp1, hp2);
  }

  // constructs segment site using either the source or the target of
  // support (that depends on the boolean is_first_exact) and the
  // intersection of support with s as the other endpoint
  Segment_Voronoi_diagram_storage_site_2(const Handle_pair& hsupport,
					 const Handle_pair& hp,
					 bool is_first_exact) {
    initialize_site(hsupport, hp, is_first_exact);
  }

#if 0
  Segment_Voronoi_diagram_storage_site_2(const Object &o) {
    if ( assign(p_, o) ) {
      initialize_site(p);
      return;
    }

    Segment_2 s;
    if ( assign(s, o) ) {
      initialize_site(s);
      return;
    }

    defined_ = false;
  }
#endif

  // PREDICATES
  //-----------
  bool is_defined() const { return type_ != 0; }
  bool is_point() const { return (type_ & 3) == 1; }
  bool is_segment() const { return (type_ & 3) == 2; }
  bool is_exact() const { return (type_ & 12) == 0; }
  bool is_exact(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    if ( i == 0 ) { return (type_ & 4) == 0; }
    return (type_ & 8) == 0;
  }

  // ACCESS METHODS
  //---------------
  Handle      point_handle() const { return h_[0]; }

  Handle_pair segment_handle() const {
    return Handle_pair(h_[0], h_[1]);
  }

  Handle_pair supporting_segment_handle() const {
    CGAL_precondition( is_segment() );
    return Handle_pair(h_[0], h_[1]);
  }

  Handle_pair supporting_segment_handle(unsigned int i) const {
    CGAL_precondition( is_point() && !is_exact() && i < 2 );
    if ( i == 0 ) {
      return Handle_pair(h_[2], h_[3]);
    } else {
      return Handle_pair(h_[4], h_[5]);
    }
  }

  Handle_pair crossing_segment_handle(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_exact() );
    CGAL_precondition( i < 2 && !is_exact(i) );
    if ( i == 0 ) {
      return Handle_pair(h_[2], h_[3]);
    } else {
      return Handle_pair(h_[4], h_[5]);
    }
  }

  // TO BE REMOVED
  void set_point(const Point_2& p) {};

  Site_2 site() const {
    if ( is_point() ) {
      if ( is_exact() ) {
	return Site_2(point());
      } else {
	return Site_2(supporting_segment(0),
		      supporting_segment(1));
      }
    } else {
      if ( is_exact() ) {
	return Site_2( segment() );
      } else if ( is_exact(0) ) {
	return Site_2( supporting_segment(), crossing_segment(1), true);
      } else if ( is_exact(1) ) {
	return Site_2( supporting_segment(), crossing_segment(0), false);
      } else {
	return Site_2( supporting_segment(), crossing_segment(0),
		       crossing_segment(1));
      }
    }
  }
protected:
  Point_2 point() const { 
    CGAL_precondition ( is_point() );
    if ( !is_exact() ) {
      return compute_intersection_point(*h_[2], *h_[3], *h_[4], *h_[5]);
    } else {
      return *h_[0];
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
      return Self(h_[0]);
    } else {
      return
	Self(Handle_pair(h_[0], h_[1]), Handle_pair(h_[2], h_[3]));
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(1) ) {
      return Self(h_[1]);
    } else {
      return
	Self(Handle_pair(h_[0], h_[1]), Handle_pair(h_[4], h_[5]));
    }
  }

  Self opposite_site() const {
    CGAL_precondition( is_segment() );

    Handle_pair supp(h_[1], h_[0]);

    if ( is_exact() ) {
      return Self(supp);
    }

    CGAL_assertion( !is_exact(0) || !is_exact(1) );

    if ( is_exact(0) && !is_exact(1) ) {
      return Self(supp, Handle_pair(h_[4], h_[5]), false);
    } else if ( !is_exact(0) && is_exact(1) ) {
      return Self(supp, Handle_pair(h_[2], h_[3]), true);
    } else {
      return Self(supp, Handle_pair(h_[4], h_[5]),
		  Handle_pair(h_[2], h_[3]));
    }
  }

  Segment_2 supporting_segment() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() ) {
      return segment();
    } else {
      return Segment_2(*h_[0], *h_[1]);
    }
  }

  Segment_2 supporting_segment(unsigned int i) const {
    //    CGAL_precondition( is_point() && !input_ && i < 2 );
    CGAL_precondition( is_point() && !is_exact() && i < 2 );
    if ( i == 0 ) {
      return Segment_2(*h_[2], *h_[3]);
    } else {
      return Segment_2(*h_[4], *h_[5]);
    }
  }

  Segment_2 crossing_segment(unsigned int i) const {
    //    CGAL_precondition( is_segment() && !input_ );
    CGAL_precondition( is_segment() && !is_exact() );
    //    CGAL_precondition( i < 2 && !is_exact_[i] );
    CGAL_precondition( i < 2 && !is_exact(i) );
    if ( i == 0 ) {
      return Segment_2(*h_[2], *h_[3]);
    } else {
      return Segment_2(*h_[4], *h_[5]);
    }
  }
public:
  // SET METHODS
  //------------
  void set_point(const Handle& h) {
    initialize_site(h);
  }

  void set_segment(const Handle_pair& hp) {
    initialize_site(hp);
  }

  void set_point(const Handle_pair& hp1, const Handle_pair& hp2) {
    initialize_site(hp1, hp2);
  }

  void set_segment(const Handle_pair& hsupport,
		   const Handle_pair& hp1,
		   const Handle_pair& hp2) {
    initialize_site(hsupport, hp1, hp2);
  }

  void set_segment(const Handle_pair& hsupport,
		   const Handle_pair& hp,
		   bool is_first_exact) {
    initialize_site(hsupport, hp, is_first_exact);
  }

public:
  std::ostream& write(std::ostream& os)
  {
    return os << (*this);
  }

protected:
  // INITIALIZATION
  //---------------
  void initialize_site(const Handle& h)
  {
    type_ = 1;
    h_[0] = h;
  }

  void initialize_site(const Handle_pair& hp)
  {
    type_ = 2;
    h_[0] = hp.first;
    h_[1] = hp.second;
  }
  void initialize_site(const Handle_pair& hp1,
		       const Handle_pair& hp2)
  {
    // MK: Sort the segments s1 and s2 in lexicographical order so
    //     that the computation of the intersection point is always
    //     done in the same manner (?)
    type_ = 5;
    h_[2] = hp1.first;
    h_[3] = hp1.second;
    h_[4] = hp2.first;
    h_[5] = hp2.second;
  }


  void initialize_site(const Handle_pair& hsupport,
		       const Handle_pair& hp1,
		       const Handle_pair& hp2)
  {
    type_ = 14;
    h_[0] = hsupport.first;
    h_[1] = hsupport.second;
    h_[2] = hp1.first;
    h_[3] = hp1.second;
    h_[4] = hp2.first;
    h_[5] = hp2.second;
  }

  void initialize_site(const Handle_pair& hsupport,
		       const Handle_pair& hp,
		       bool is_first_exact)
  {
    type_ = (is_first_exact ? 10 : 6);
    h_[0] = hsupport.first;
    h_[1] = hsupport.second;
    if ( is_first_exact ) {
      h_[4] = hp.first;
      h_[5] = hp.second;
    } else {
      h_[2] = hp.first;
      h_[3] = hp.second;
    }
  }

  // CONSTRUCTION METHODS
  //---------------------
  Point_2 compute_source() const {
    CGAL_precondition( is_segment() );
    //    if ( input_ || is_exact_[0] ) {
    if ( is_exact() || is_exact(0) ) {
      return *h_[0];
    } else {
      return compute_intersection_point(*h_[0], *h_[1], *h_[2], *h_[3]);
    }
  }

  Point_2 compute_target() const {
    CGAL_precondition( is_segment() );
    //    if ( input_ || is_exact_[1] ) {
    if ( is_exact() || is_exact(1) ) {
      return *h_[1];
    } else {
      return compute_intersection_point(*h_[0], *h_[1], *h_[4], *h_[5]);
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
  Handle h_[6];
  char type_;
};

//-------------------------------------------------------------------------

template <class R, class H>
std::ostream&
operator<<(std::ostream& os, 
	   const Segment_Voronoi_diagram_storage_site_2<R,H>& s)
{
  if (!s.is_defined())
    return os << "u";
  if (s.is_point())
    return os << "p " << s.point ();
  return os << "s " << s.segment ();
}

template < class R, class H, class Stream >
Stream&
operator<<(Stream& str, Segment_Voronoi_diagram_storage_site_2<R,H>& t)
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

#endif // SEGMENT_VORONOI_DIAGRAM_STORAGE_SITE_H
