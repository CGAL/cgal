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



#ifndef SEGMENT_VORONOI_DIAGRAM_SIMPLE_STORAGE_SITE_H
#define SEGMENT_VORONOI_DIAGRAM_SIMPLE_STORAGE_SITE_H

#include <iostream>
#include <CGAL/assertions.h>


CGAL_BEGIN_NAMESPACE

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class Gt, class H>
class Segment_Voronoi_diagram_simple_storage_site_2 
{
public:
  typedef Gt                             Geom_traits;
  typedef H                              Point_handle;
  typedef typename Geom_traits::Site_2   Site_2;

protected:
  typedef Point_handle                   Handle;

  typedef
  Segment_Voronoi_diagram_simple_storage_site_2<Geom_traits,Handle>
  Self;

public:
  Segment_Voronoi_diagram_simple_storage_site_2() : type_(0) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_simple_storage_site_2(const Handle& hp) {
    initialize_site(hp);
  }

  // constructs segment site using input segment
  Segment_Voronoi_diagram_simple_storage_site_2(const Handle& hp1,
						const Handle& hp2) {
    initialize_site(hp1, hp2);
  }

  // the compiler complains that it cannot find this constructor;
  // solution: make the insert_intersecting_segment a template
  // method...
  template<class A1, class A2, class A3, class A4>
  Segment_Voronoi_diagram_simple_storage_site_2(const A1&, const A2&,
						const A3&, const A4&) {
    CGAL_assertion( false );
  }

  template<class A1, class A2, class A3, class A4, class A5>
  Segment_Voronoi_diagram_simple_storage_site_2(const A1&, const A2&,
						const A3&, const A4&,
						const A5&) {
    CGAL_assertion( false );
  }

  template<class A1, class A2, class A3, class A4, class A5, class A6>
  Segment_Voronoi_diagram_simple_storage_site_2(const A1&, const A2&,
						const A3&, const A4&,
						const A5&, const A6&) {
    CGAL_assertion( false );
  }

public:
  // PREDICATES
  //-----------
  bool is_defined() const { return type_; }
  bool is_point() const { return type_ == 1; }
  bool is_segment() const { return type_ == 2; }
  bool is_exact() const { return true; }
  bool is_exact(unsigned int i) const { return true; }

  // ACCESS METHODS
  //---------------
  const Handle& point_handle(unsigned int i) const {
    CGAL_precondition( i < 6 );
    return h_[i];
  }

  Self supporting_segment_site() const {
    CGAL_precondition( is_segment() );
    return Self(h_[0], h_[1]);
  }

  Self supporting_segment_site(unsigned int i) const {
    CGAL_precondition( false );
    CGAL_precondition( is_point() && i < 2 );
    return Self(h_[0], h_[0]);
  }

  Self crossing_segment_handle(unsigned int i) const {
    CGAL_precondition( false );
    CGAL_precondition( is_segment() && i < 2 );
    return Self(h_[0], h_[1]);
  }

  Site_2 site() const {
    if ( is_point() ) {
      return Site_2(*h_[0]);
    } else {
      return Site_2(*h_[0], *h_[1]);
    }
  }

public:
  // SET METHODS
  //------------
  void set_point(const Handle& hp) {
    initialize_site(hp);
  }

  void set_segment(const Handle& hp1, const Handle& hp2) {
    initialize_site(hp1, hp2);
  }

  // the compiler complains that it cannot find this constructor;
  // solution: make the insert_intersecting_segment a template
  // method...
  template<class A1, class A2, class A3, class A4>
  void set_point(const A1&, const A2&, const A3&, const A4&) {
    CGAL_assertion(false);
  }

  template<class A1, class A2, class A3, class A4, class A5>
  void set_segment(const A1&, const A2&, const A3&, const A4&,
		   const A5&) {
    CGAL_assertion(false);
  }

  template<class A1, class A2, class A3, class A4, class A5, class A6>
  void set_segment(const A1&, const A2&, const A3&, const A4&,
		   const A5&, const A6&) {
    CGAL_assertion(false);
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

protected:
  Handle h_[2];
  char type_;
};

//-------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_SIMPLE_STORAGE_SITE_H
