// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_2_H 1

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <set>
#include <CGAL/Segment_Delaunay_graph_storage_site_2.h>
#include <CGAL/Segment_Delaunay_graph_simple_storage_site_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Construct_storage_site_2.h>


namespace CGAL {

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace SegmentDelaunayGraph_2 {

namespace Internal {

  template<class Gt, class USE_SIMPLE_STORAGE_SITE_Tag>
  struct Which_storage_site;

  // use the simple storage site
  template<class Gt>
  struct Which_storage_site<Gt,Tag_false>
  {
    typedef Gt         Geom_traits;
    typedef Tag_false  Storage_site_tag;

    typedef Segment_Delaunay_graph_simple_storage_site_2<Geom_traits>
    Storage_site_2;
  };

  // use the full storage site
  template<class Gt>
  struct Which_storage_site<Gt,Tag_true>
  {
    typedef Gt         Geom_traits;
    typedef Tag_true   Storage_site_tag;

    typedef Segment_Delaunay_graph_storage_site_2<Geom_traits>
    Storage_site_2;
  };


} // namespace Internal

} //namespace SegmentDelaunayGraph_2

//----------------------------------------------------------------------

template<class Gt>
class Segment_Delaunay_graph_storage_traits_2
{
public:
  typedef Gt                                       Geom_traits;
  typedef typename Geom_traits::Point_2            Point_2;
  typedef typename Geom_traits::Site_2             Site_2;
  typedef std::set<Point_2>                        Point_container;
  typedef typename Point_container::iterator       Point_handle;
  typedef typename Point_container::const_iterator const_Point_handle;

private:
  typedef Segment_Delaunay_graph_storage_traits_2<Geom_traits>   Self;
  typedef typename Geom_traits::Intersections_tag                ITag;

public:
  typedef typename
  CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Internal::
  Which_storage_site<Self,ITag>::Storage_site_2
  Storage_site_2;

  typedef
  CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Construct_storage_site_2<Self>
  Construct_storage_site_2;

  // MK::FIGURE OUT HOW TO PASS A REFERENCE TO GEOM_TRAITS AND HAVE
  // DEFAULT CONSTRUCTOR AS WELL IF POSSIBLE
  Segment_Delaunay_graph_storage_traits_2
  (const Geom_traits& gt = Geom_traits()) : gt_(gt) {}

  inline const Geom_traits& geom_traits() const { return gt_; }

  inline Construct_storage_site_2
  construct_storage_site_2_object() const {
    return Construct_storage_site_2();
  }

private:
  Geom_traits gt_;
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_STORAGE_TRAITS_2_H
