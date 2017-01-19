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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCT_STORAGE_SITE_WITH_INFO_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCT_STORAGE_SITE_WITH_INFO_2_H 1

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <CGAL/Segment_Delaunay_graph_2/Construct_storage_site_2.h>


namespace CGAL {

namespace SegmentDelaunayGraph_2 {

template<class STraits>
class Construct_storage_site_with_info_2
  : public Construct_storage_site_2<STraits>
{
public:
  typedef STraits                                    Storage_traits;
  typedef typename Storage_traits::Storage_site_2    Storage_site_2;
  typedef typename Storage_traits::Point_handle      Point_handle;

  typedef Storage_site_2                             result_type;

protected:
  typedef Construct_storage_site_2<Storage_traits>   Base;
  typedef typename Storage_traits::Info              Info;
  typedef typename Storage_traits::Convert_info      Convert_info;
  typedef typename Storage_traits::Merge_info        Merge_info;

public:
  // constructs the point of intersection
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1) const {
    Storage_site_2 ssx = Base::operator()(ss0, ss1);
    Info infox = Merge_info()(ss0.info(), ss1.info());
    ssx.set_info(infox);
    return ssx;
  }

  // constructs the subsegment with supporting segment ss0 and
  // endpoints the point of intersection of ss1 and ss0; the boolean
  // determines if the first or segment subsegment is constructed
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1,
			 bool first) const {
    Storage_site_2 s = Base::operator()(ss0, ss1, first);
    Info is = Convert_info()(ss0.info(), ss1.info(), first);
    s.set_info(is);
    return s;
  }

  using Base::operator();
};


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCT_STORAGE_SITE_WITH_INFO_2_H
