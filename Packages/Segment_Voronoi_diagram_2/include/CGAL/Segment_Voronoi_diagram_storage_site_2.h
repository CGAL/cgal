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
  //  typedef typename R::Point_2   Point_2;
  //  typedef typename R::Segment_2 Segment_2;
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

public:
  // PREDICATES
  //-----------
  bool is_defined() const { return type_; }
  bool is_point() const { return (type_ & 3) == 1; }
  bool is_segment() const { return (type_ & 3) == 2; }
  bool is_exact() const { return !(type_ & 12); }
  bool is_exact(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    if ( i == 0 ) { return !(type_ & 4); }
    return !(type_ & 8);
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

  Site_2 site() const {
    if ( is_point() ) {
      if ( is_exact() ) {
	return Site_2(*h_[0]);
      } else {
	return Site_2(*h_[2], *h_[3], *h_[4], *h_[5]);
      }
    } else {
      if ( is_exact() ) {
	return Site_2(*h_[0], *h_[1]);
      } else if ( is_exact(0) ) {
	return Site_2(*h_[0], *h_[1], *h_[4], *h_[5], true);
      } else if ( is_exact(1) ) {
	return Site_2(*h_[0], *h_[1], *h_[2], *h_[3], false);
      } else {
	return Site_2(*h_[0], *h_[1], *h_[2], *h_[3], *h_[4], *h_[5]);
      }
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

protected:
  Handle h_[6];
  char type_;
};

//-------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_STORAGE_SITE_H
