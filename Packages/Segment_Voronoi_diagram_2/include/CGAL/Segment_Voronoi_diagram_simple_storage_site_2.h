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

template <class R_, class H_>
class Segment_Voronoi_diagram_simple_storage_site_2 
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
  typedef Segment_Voronoi_diagram_simple_storage_site_2<Rep,Handle>  Self;

public:
  Segment_Voronoi_diagram_simple_storage_site_2() : type_(0) {}

  // constructs point site using input point
  Segment_Voronoi_diagram_simple_storage_site_2(const Handle &h) {
    initialize_site(h);
  }

  // constructs segment site using input segment
  Segment_Voronoi_diagram_simple_storage_site_2(const Handle_pair &hp) {
    initialize_site(hp);
  }

  // the compiler complains that it cannot find this constructor;
  // solution: make the insert_intersecting_segment a template
  // method...
  template<class A1, class A2>
  Segment_Voronoi_diagram_simple_storage_site_2(const A1&, const A2&) {
    CGAL_assertion( false );
  }

  template<class A1, class A2, class A3>
  Segment_Voronoi_diagram_simple_storage_site_2(const A1&, const A2&,
						const A3&) {
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
  Handle      point_handle() const { return h_[0]; }

  Handle_pair segment_handle() const {
    return Handle_pair(h_[0], h_[1]);
  }

  Handle_pair supporting_segment_handle() const {
    CGAL_precondition( is_segment() );
    return Handle_pair(h_[0], h_[1]);
  }

  Handle_pair supporting_segment_handle(unsigned int i) const {
    CGAL_precondition( is_point() && i < 2 );
    return Handle_pair(h_[0], h_[0]);
  }

  Handle_pair crossing_segment_handle(unsigned int i) const {
    CGAL_precondition( is_segment() && i < 2 );
    return Handle_pair(h_[0], h_[1]);
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
  void set_point(const Handle& h) {
    initialize_site(h);
  }

  void set_segment(const Handle_pair& hp) {
    initialize_site(hp);
  }

  // the compiler complains that it cannot find this constructor;
  // solution: make the insert_intersecting_segment a template
  // method...
  template<class A1, class A2>
  void set_point(const A1&, const A2&)
  {
    CGAL_assertion(false);
  }

  template<class A1, class A2, class A3>
  void set_segment(const A1&, const A2&, const A3&)
  {
    CGAL_assertion(false);
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

protected:
  Handle h_[2];
  char type_;
};

//-------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_SIMPLE_STORAGE_SITE_H
