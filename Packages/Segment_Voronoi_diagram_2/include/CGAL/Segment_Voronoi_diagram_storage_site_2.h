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

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

CGAL_BEGIN_NAMESPACE

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class Gt>
class Segment_Voronoi_diagram_storage_site_2 
{
public:
  typedef Gt                             Geom_traits;
  typedef typename Geom_traits::Site_2   Site_2;
  typedef typename std::list<typename Site_2::Point_2>::iterator Point_handle;

protected:
  typedef Point_handle                   Handle;

  typedef Segment_Voronoi_diagram_storage_site_2<Geom_traits>  Self;

public:
  Segment_Voronoi_diagram_storage_site_2() : type_(0) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_storage_site_2(const Handle& hp) {
    initialize_site(hp);
  }

  // constructs segment site corresponding to the segment (*hp1,*hp2)
  Segment_Voronoi_diagram_storage_site_2(const Handle& hp1,
					 const Handle& hp2) {
    initialize_site(hp1, hp2);
  }

  // constructs point site using the point of intersection of the
  // segments (*hp1,*hp2) and (*hq1,*hq2)
  Segment_Voronoi_diagram_storage_site_2(const Handle& hp1,
					 const Handle& hp2,
					 const Handle& hq1,
					 const Handle& hq2) {
    initialize_site(hp1, hp2, hq1, hq2);
  }

  // constructs segment site whose endpoints are the points of
  // intersection of the pairs of segments (*hp1,*hp2), (*hq1,*hq2)
  // and (*hp1,*hp2), (*hr1,*hr2)
  Segment_Voronoi_diagram_storage_site_2(const Handle& hp1,
					 const Handle& hp2,
					 const Handle& hq1,
					 const Handle& hq2,
					 const Handle& hr1,
					 const Handle& hr2) {
    initialize_site(hp1, hp2, hq1, hq2, hr1, hr2);
  }

  // constructs segment site using either the source or the target of
  // (*hp1,*hp2) (that depends on the boolean is_first_exact) and the
  // intersection of (*hp1,*hp2) with (*hq1,*hq2) as the other endpoint
  Segment_Voronoi_diagram_storage_site_2(const Handle& hp1,
					 const Handle& hp2,
					 const Handle& hq1,
					 const Handle& hq2,
					 bool is_first_exact) {
    initialize_site(hp1, hp2, hq1, hq2, is_first_exact);
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
  const Handle& point_handle(unsigned int i) const {
    CGAL_precondition( i < 6 );
    return h_[i];
  }

  Self source_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(0) ) {
      return Self(h_[0]);
    } else {
      return Self(h_[0], h_[1], h_[2], h_[3]);
    }
  }

  Self target_site() const {
    CGAL_precondition( is_segment() );
    if ( is_exact() || is_exact(1) ) {
      return Self(h_[1]);
    } else {
      return Self(h_[0], h_[1], h_[4], h_[5]);
    }
  }

  Self supporting_segment_site() const {
    CGAL_precondition( is_segment() );
    return Self(h_[0], h_[1]);
  }

  Self supporting_segment_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_exact() && i < 2 );
    if ( i == 0 ) {
      return Self(h_[2], h_[3]);
    } else {
      return Self(h_[4], h_[5]);
    }
  }

  Self crossing_segment_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_exact() );
    CGAL_precondition( i < 2 && !is_exact(i) );
    if ( i == 0 ) {
      return Self(h_[2], h_[3]);
    } else {
      return Self(h_[4], h_[5]);
    }
  }

  Site_2 site() const {
    if ( is_point() ) {
      if ( is_exact() ) {
	return Site_2::construct_site_2(*h_[0]);
      } else {
	return Site_2::construct_site_2(*h_[2], *h_[3], *h_[4], *h_[5]);
      }
    } else {
      if ( is_exact() ) {
	return Site_2::construct_site_2(*h_[0], *h_[1]);
      } else if ( is_exact(0) ) {
	return Site_2::construct_site_2(*h_[0], *h_[1], *h_[4],
					*h_[5], true);
      } else if ( is_exact(1) ) {
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

protected:
  Handle h_[6];
  char type_;
};

//-------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_STORAGE_SITE_H
