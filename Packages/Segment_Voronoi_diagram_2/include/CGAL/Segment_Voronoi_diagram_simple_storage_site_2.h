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

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

CGAL_BEGIN_NAMESPACE

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

template <class Gt>
class Segment_Voronoi_diagram_simple_storage_site_2 
{
public:
  typedef Gt                             Geom_traits;
  typedef typename Geom_traits::Site_2   Site_2;
  typedef typename std::set<typename Site_2::Point_2>::iterator Point_handle;
  //  typedef typename std::list<typename Site_2::Point_2>::iterator Point_handle;

protected:
  typedef Point_handle                   Handle;

  typedef
  Segment_Voronoi_diagram_simple_storage_site_2<Geom_traits>
  Self;

public:
  // constructs point site using input point
  static Self construct_storage_site_2(const Handle& hp) {
    Self t;
    t.initialize_site(hp);
    return t;
  }

  // constructs segment site using input segment
  static Self construct_storage_site_2(const Handle& hp1,
				       const Handle& hp2) {
    Self t;
    t.initialize_site(hp1, hp2);
    return t;
  }

public:
  // DEFAULT CONSTRUCTOR
  //--------------------
  Segment_Voronoi_diagram_simple_storage_site_2() : type_(0) {}

  // COPY CONSTRUCTOR
  //-----------------
  Segment_Voronoi_diagram_simple_storage_site_2(const Self& other) {
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
  bool is_defined() const { return type_; }
  bool is_point() const { return type_ == 1; }
  bool is_segment() const { return type_ == 2; }
  bool is_input() const { return true; }
  bool is_input(unsigned int i) const { return true; }

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

  // the following methods should never be called; they have been
  // defined in order for this class to be a model of the
  // SegmentVoronoiDiagramStorageSite_2 concept.
  const Handle& source_of_supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() );
    return h_[0];
  }

  const Handle& target_of_supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() );
    return h_[0];
  }

  const Handle& source_of_crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input(i) );
    return h_[0];
  }

  const Handle& target_of_crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input(i) );
    return h_[0];
  }

  Site_2 site() const {
    if ( is_point() ) {
      return Site_2::construct_site_2(*h_[0]);
    } else {
      return Site_2::construct_site_2(*h_[0],*h_[1]);
    }
  }

  Self source_site() const
  {
    CGAL_precondition( is_segment() );
    return Self(h_[0]);
  }

  Self target_site() const
  {
    CGAL_precondition( is_segment() );
    return Self(h_[1]);
  }

  Self supporting_site() const {
    CGAL_precondition( is_segment() );
    return Self(h_[0], h_[1]);
  }

  Self supporting_site(unsigned int i) const {
    CGAL_precondition( is_point() && !is_input() && i < 2 );
    return Self(h_[0], h_[0]);
  }

  Self crossing_site(unsigned int i) const {
    CGAL_precondition( is_segment() && !is_input() );
    CGAL_precondition( i < 2 && !is_input(i) );
    return Self(h_[0], h_[0]);
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

  void copy_from(const Self& other) {
    type_ = other.type_;

    if ( !other.is_defined() ) { return; }

    h_[0] = other.h_[0];

    if ( other.is_segment() ) {
      h_[1] = other.h_[1];
    }
  }

protected:
  Handle h_[2];
  char type_;
};

//-------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_SIMPLE_STORAGE_SITE_H
