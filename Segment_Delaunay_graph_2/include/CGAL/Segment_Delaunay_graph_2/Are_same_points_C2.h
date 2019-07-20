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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARE_SAME_POINTS_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARE_SAME_POINTS_C2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

template<class K>
class Are_same_points_C2
{
private:
  typedef typename K::Point_2     Point_2;
  typedef typename K::Segment_2   Segment_2;
  typedef typename K::Site_2      Site_2;
  typedef typename K::Compare_x_2 Compare_x_2;
  typedef typename K::Compare_y_2 Compare_y_2;
  typedef typename K::Boolean     Boolean;

  typedef typename K::Intersections_tag  ITag;

  Compare_x_2 compare_x_2;
  Compare_y_2 compare_y_2;

  Boolean   are_same(const Point_2& p, const Point_2& q) const
  {
    return
      compare_x_2(p, q) == EQUAL && compare_y_2(p, q) == EQUAL;
  }

  Boolean   are_same(const Site_2& s, const Site_2& t) const
  {
    return
      ( are_same(s.source(), t.source()) &&
	are_same(s.target(), t.target()) ) ||
      ( are_same(s.source(), t.target()) &&
	are_same(s.target(), t.source()) );
  }

  Boolean   predicate(const Site_2& p, const Site_2& q, const Tag_false&) const
  {
    return are_same(p.point(), q.point()); 
  }

  Boolean   predicate(const Site_2& p, const Site_2& q, const Tag_true&) const
  {
    if ( !p.is_input() && !q.is_input() ) {
      Site_2 s[2] = { p.supporting_site(0), p.supporting_site(1) };
      Site_2 t[2] = { q.supporting_site(0), q.supporting_site(1) };

      if (  ( are_same(s[0], t[0]) && are_same(s[1], t[1]) ) ||
	    ( are_same(s[0], t[1]) && are_same(s[1], t[0]) )  ) {
	return true;
      }
    }

    return predicate(p, q, Tag_false());
  }

public:
  typedef Boolean        result_type;
  typedef Site_2         argument_type;

  Boolean   operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    return predicate(p, q, ITag());
  }
};

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARE_SAME_POINTS_C2_H
