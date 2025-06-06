// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_COMPARE_Y_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_COMPARE_Y_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {


//-----------------------------------------------------------------------
//                           compare y
//-----------------------------------------------------------------------

template< class K >
class Compare_y_2
{
public:
  typedef typename K::Site_2                Site_2;
  typedef typename K::Point_2               Point_2;
  typedef typename K::Comparison_result     Comparison_result;

private:
  typedef typename K::Compare_y_2           Kernel_compare_y_2;

public:
  inline
  Comparison_result operator()(const Point_2& p, const Point_2& q) const
  {
    return Kernel_compare_y_2()( p, q );
  }

  inline
  Comparison_result operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    return Kernel_compare_y_2()( p.point(), q.point() );
  }
};

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_COMPARE_Y_2_H
