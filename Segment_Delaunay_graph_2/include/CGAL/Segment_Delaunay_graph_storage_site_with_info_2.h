// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_SITE_WITH_INFO_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_SITE_WITH_INFO_2_H 1

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <iostream>
#include <CGAL/assertions.h>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/Segment_Delaunay_graph_2/Construct_storage_site_2.h>

namespace CGAL {

  /** A Site is either a point or a segment or a point defined as the
      intersection of two non-parallel segments (if defined)
   */

namespace SegmentDelaunayGraph_2 {

template<class STraits> class Construct_storage_site_with_info_2;

} //namespace SegmentDelaunayGraph_2


template <class STraits, typename Info_, class Base_storage_site>
class Segment_Delaunay_graph_storage_site_with_info_2
  : public Base_storage_site
{
  typedef Base_storage_site                       Base;

  friend class
  CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Construct_storage_site_2<STraits>;

public:
  typedef STraits                                 Storage_traits;
  typedef Info_                                   Info;
  typedef typename Storage_traits::Geom_traits    Geom_traits;
  typedef typename Geom_traits::Site_2            Site_2;
  typedef typename Storage_traits::Point_handle   Point_handle;

  struct Has_info_tag {};

protected:
  typedef Point_handle                            Handle;

  typedef
  Segment_Delaunay_graph_storage_site_with_info_2<Storage_traits,Info,Base>
  Self;

public:
  // DEFAULT CONSTRUCTOR
  //--------------------
  Segment_Delaunay_graph_storage_site_with_info_2()
    : Base(), info_set_(false) {}

private:
  // CONSTRUCTOR FROM BASE SITE
  //---------------------------
  Segment_Delaunay_graph_storage_site_with_info_2(const Base& base)
    : Base(base), info_set_(false) {}

public:
  // COPY CONSTRUCTOR
  //-----------------
  Segment_Delaunay_graph_storage_site_with_info_2(const Self& other) 
    : Base(other)
  {
    info_set_ = other.info_set_;
    if ( info_set_ ) {
      info_ = other.info_;
    }
  }

  // ASSIGNMENT OPERATOR
  //--------------------
  Self& operator=(const Self& other) {
    Base::operator=(other);
    info_set_ = other.info_set_;
    if ( info_set_ ) {
      info_ = other.info_;
    }
    return *this;
  }

public:
  inline Self source_site() const {
    return Base::source_site();
  }

  inline Self target_site() const {
    return Base::target_site();
  }

  inline Self supporting_site() const {
    return Base::supporting_site();
  }

  inline Self supporting_site(unsigned int i) const {
    return Base::supporting_site(i);
  }

  inline Self crossing_site(unsigned int i) const {
    return Base::crossing_site(i);
  }

  // ACCESS TO INFO
  //---------------
  inline const Info& info() const {
    CGAL_precondition( info_has_been_set() );
    return info_;
  }

  inline void set_info(const Info& info) {
    info_set_ = true;
    info_ = info;
  }

  inline bool info_has_been_set() const { return info_set_; }

protected:
  bool info_set_;
  Info info_;
};

//-------------------------------------------------------------------------

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_SITE_WITH_INFO_2_H
