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

  Site_2 site() const {
    if ( is_point() ) {
#ifdef USE_SC
      return Site_2(*h_[0]);
#else
      return Site_2::construct_site_2(*h_[0]);
#endif
    } else {
#ifdef USE_SC
      return Site_2(*h_[0], *h_[1]);
#else
      return Site_2::construct_site_2(*h_[0],*h_[1]);
#endif
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

protected:
  Handle h_[2];
  char type_;
};

//-------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // SEGMENT_VORONOI_DIAGRAM_SIMPLE_STORAGE_SITE_H
