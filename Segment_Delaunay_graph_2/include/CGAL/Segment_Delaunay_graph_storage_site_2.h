// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_SITE_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_SITE_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <iostream>
#include <CGAL/assertions.h>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

namespace SegmentDelaunayGraph_2 {

template<class STraits> class Construct_storage_site_2;

} //namespace SegmentDelaunayGraph_2


template <class STraits>
class Segment_Delaunay_graph_storage_site_2 
{
  friend class
  CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Construct_storage_site_2<STraits>;

public:
  typedef STraits                                 Storage_traits;
  typedef typename Storage_traits::Geom_traits    Geom_traits;
  typedef typename Geom_traits::Site_2            Site_2;
  typedef typename Storage_traits::Point_handle   Point_handle;

protected:
  typedef Point_handle                   Handle;

  typedef Segment_Delaunay_graph_storage_site_2<Storage_traits>  Self;

protected:
  // constructs point site using input point
  static Self construct_storage_site_2(const Handle& hp) {
    Self t;
    t.initialize_site(hp);
    return t;
  }

  // constructs segment site corresponding to the segment (*hp1,*hp2)
  static Self construct_storage_site_2(const Handle& hp1,
				       const Handle& hp2) {
    Self t;
    t.initialize_site(hp1, hp2);
    return t;
  }

  // constructs point site using the point of intersection of the
  // segments (*hp1,*hp2) and (*hq1,*hq2)
  static Self construct_storage_site_2(const Handle& hp1,
				       const Handle& hp2,
				       const Handle& hq1,
				       const Handle& hq2) {
    Self t;
    t.initialize_site(hp1, hp2, hq1, hq2);
    return t;
  }

  // constructs segment site whose endpoints are the points of
  // intersection of the pairs of segments (*hp1,*hp2), (*hq1,*hq2)
  // and (*hp1,*hp2), (*hr1,*hr2)
  static Self construct_storage_site_2(const Handle& hp1,
				       const Handle& hp2,
				       const Handle& hq1,
				       const Handle& hq2,
				       const Handle& hr1,
				       const Handle& hr2) {
    Self t;
    t.initialize_site(hp1, hp2, hq1, hq2, hr1, hr2);
    return t;
  }

  // constructs segment site using either the source or the target of
  // (*hp1,*hp2) (that depends on the boolean is_first_exact) and the
  // intersection of (*hp1,*hp2) with (*hq1,*hq2) as the other endpoint
  static Self construct_storage_site_2(const Handle& hp1,
				       const Handle& hp2,
				       const Handle& hq1,
				       const Handle& hq2,
				       bool is_first_exact) {
    Self t;
    t.initialize_site(hp1, hp2, hq1, hq2, is_first_exact);
    return t;
  }


public:
  // DEFAULT CONSTRUCTOR
  //--------------------
  Segment_Delaunay_graph_storage_site_2() : type_(0) {}

  // COPY CONSTRUCTOR
  //-----------------
  Segment_Delaunay_graph_storage_site_2(const Self& other) {
    copy_from(other);
  }

  // ASSIGNMENT OPERATOR
  //--------------------
  Self& operator=(const Self& other) {
    copy_from(other);
    return *this;
  }

public:
  // PREDICATES
  //-----------
  bool is_defined() const { return type_ != 0; }
  bool is_point() const { return (type_ & 3) == 1; }
  bool is_segment() const { return (type_ & 3) == 2; }
  bool is_input() const { return !(type_ & 12); }
  bool is_input(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    if ( i == 0 ) { return !(type_ & 4); }
    return !(type_ & 8);
  }

  // ACCESS METHODS
  //---------------
  const Handle& point() const {
    CGAL_precondition( is_point() && is_input() );
    return h_[0];
  }

  const Handle& source_of_supporting_site() const {
    CGAL_precondition( is_segment() );
    return h_[0];
  }

  const Handle& target_of_supporting_site() const {
    CGAL_precondition( is_segment() );
    return h_[1];
  }

  const Handle& source_of_supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() );
    return (i == 0) ? h_[2] : h_[4];
  }

  const Handle& target_of_supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() );
    return (i == 0) ? h_[3] : h_[5];
  }

  const Handle& source_of_crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input(i) );
    return (i == 0) ? h_[2] : h_[4];
  }

  const Handle& target_of_crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input(i) );
    return (i == 0) ? h_[3] : h_[5];
  }

  Self source_site() const {
    CGAL_precondition( is_segment() );
    if ( is_input() || is_input(0) ) {
      return construct_storage_site_2(h_[0]);
    } else {
      return construct_storage_site_2(h_[0], h_[1], h_[2], h_[3]);
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( is_input() || is_input(1) ) {
      return construct_storage_site_2(h_[1]);
    } else {
      return construct_storage_site_2(h_[0], h_[1], h_[4], h_[5]);
    }
  }

  Self supporting_site() const {
    CGAL_precondition( is_segment() );
    return construct_storage_site_2(h_[0], h_[1]);
  }

  Self supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() && i < 2 );
    if ( i == 0 ) {
      return construct_storage_site_2(h_[2], h_[3]);
    } else {
      return construct_storage_site_2(h_[4], h_[5]);
    }
  }

  Self crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input() );
    CGAL_precondition( i < 2 && !is_input(i) );
    if ( i == 0 ) {
      return construct_storage_site_2(h_[2], h_[3]);
    } else {
      return construct_storage_site_2(h_[4], h_[5]);
    }
  }

  Site_2 site() const {
    if ( is_point() ) {
      if ( is_input() ) {
	return Site_2::construct_site_2(*h_[0]);
      } else {
	return Site_2::construct_site_2(*h_[2], *h_[3], *h_[4], *h_[5]);
      }
    } else {
      if ( is_input() ) {
	return Site_2::construct_site_2(*h_[0], *h_[1]);
      } else if ( is_input(0) ) {
	return Site_2::construct_site_2(*h_[0], *h_[1], *h_[4],
					*h_[5], true);
      } else if ( is_input(1) ) {
	return Site_2::construct_site_2(*h_[0], *h_[1], *h_[2],
					*h_[3], false);
      } else {
	return Site_2::construct_site_2(*h_[0], *h_[1], *h_[2],
					*h_[3], *h_[4], *h_[5]);
      }
    }
  }

protected:
  // INITIALIZATION
  //---------------
  void initialize_site(const Handle& hp)
  {
    type_ = 1;
    h_[0] = hp;
  }

  void initialize_site(const Handle& hp1, const Handle& hp2)
  {
    type_ = 2;
    h_[0] = hp1;
    h_[1] = hp2;
  }
  void initialize_site(const Handle& hp1, const Handle& hp2,
		       const Handle& hq1, const Handle& hq2)
  {
    // MK: Sort the segments s1 and s2 in lexicographical order so
    //     that the computation of the intersection point is always
    //     done in the same manner (?)
    type_ = 5;
    h_[2] = hp1;
    h_[3] = hp2;
    h_[4] = hq1;
    h_[5] = hq2;
  }


  void initialize_site(const Handle& hp1, const Handle& hp2,
		       const Handle& hq1, const Handle& hq2,
		       const Handle& hr1, const Handle& hr2)
  {
    type_ = 14;
    h_[0] = hp1;
    h_[1] = hp2;
    h_[2] = hq1;
    h_[3] = hq2;
    h_[4] = hr1;
    h_[5] = hr2;
  }

  void initialize_site(const Handle& hp1, const Handle& hp2,
		       const Handle& hq1, const Handle& hq2,
		       bool is_first_exact)
  {
    type_ = (is_first_exact ? 10 : 6);
    h_[0] = hp1;
    h_[1] = hp2;
    if ( is_first_exact ) {
      h_[4] = hq1;
      h_[5] = hq2;
    } else {
      h_[2] = hq1;
      h_[3] = hq2;
    }
  }

  void copy_from(const Self& other) {
    type_ = other.type_;

    if ( !other.is_defined() ) { return; }

    if ( other.is_point() ) {
      if ( other.is_input() ) {
	h_[0] = other.h_[0];
      } else {
	h_[2] = other.h_[2];
	h_[3] = other.h_[3];
	h_[4] = other.h_[4];
	h_[5] = other.h_[5];
      }
    } else {
      h_[0] = other.h_[0];
      h_[1] = other.h_[1];
      if ( !other.is_input() ) {
	if ( other.is_input(0) ) {
	  h_[4] = other.h_[4];
	  h_[5] = other.h_[5];
	} else if ( other.is_input(1) ) {
	  h_[2] = other.h_[2];
	  h_[3] = other.h_[3];
	} else {
	  h_[2] = other.h_[2];
	  h_[3] = other.h_[3];
	  h_[4] = other.h_[4];
	  h_[5] = other.h_[5];
	}
      }
    }
  }

protected:
  Handle h_[6];
  char type_;
};

//-------------------------------------------------------------------------

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_SITE_H
