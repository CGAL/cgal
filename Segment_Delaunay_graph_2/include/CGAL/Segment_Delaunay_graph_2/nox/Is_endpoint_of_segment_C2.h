// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France) and
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_IS_ENDPOINT_OF_SEGMENT_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_IS_ENDPOINT_OF_SEGMENT_C2_H

#include <CGAL/Segment_Delaunay_graph_2/nox/Are_same_points_C2.h>

CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

template<class K>
class Is_endpoint_of_segment_C2
{
private:
  typedef Are_same_points_C2<K>       Are_same_points_2;

private:
  Are_same_points_2  same_points;

public:
  typedef typename K::Site_2    Site_2;
  typedef typename K::Point_2   Point_2;
  typedef bool                  result_type;

  // check whether the point p is an endpoint of the segment s
  inline
  bool operator()(const Point_2& p, const Site_2& s) const
  {
    CGAL_precondition( s.is_segment() );
#ifdef CGAL_SDG_SORT_POINTS_IN_SITE2
    const Point_2& s1 = s.source();
    Comparison_result rs_x = CGAL::compare( p.x(), s1.x() );

    if ( rs_x == SMALLER ) return false;
    if ( rs_x == EQUAL ) {
      Comparison_result rs_y = CGAL::compare(p.y(), s1.y());
      if ( rs_y == SMALLER ) { return false; }
      if ( rs_y == EQUAL ) { return true; }
    }

    return same_points(p, s.target());
#else
    return ( same_points(p, s.source()) ||
	     same_points(p, s.target()) );
#endif
  }

  inline
  bool operator()(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );

    return operator()(p.point(), s);
  }
};

CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_IS_ENDPOINT_OF_SEGMENT_C2_H
