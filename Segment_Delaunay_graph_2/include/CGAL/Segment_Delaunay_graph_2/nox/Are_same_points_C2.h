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


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_ARE_SAME_POINTS_C2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_ARE_SAME_POINTS_C2_H

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_SEGMENT_DELAUNAY_GRAPH_2_BEGIN_NAMESPACE

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

public:
  typedef Boolean        result_type;
  typedef Site_2         argument_type;

  inline
  Boolean operator()(const Point_2& p, const Point_2& q) const
  {
    return
      compare_x_2(p, q) == EQUAL && compare_y_2(p, q) == EQUAL;
  }

  inline
  Boolean operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );

    return operator()(p.point(), q.point());
  }
};

CGAL_SEGMENT_DELAUNAY_GRAPH_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_ARE_SAME_POINTS_C2_H
