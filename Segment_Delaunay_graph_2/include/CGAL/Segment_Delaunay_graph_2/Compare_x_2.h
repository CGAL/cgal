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


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_COMPARE_X_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_COMPARE_X_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

namespace CGAL {

namespace SegmentDelaunayGraph_2 {

//-----------------------------------------------------------------------
//                           compare x
//-----------------------------------------------------------------------

template< class K >
class Compare_x_2
{
public:
  typedef typename K::Site_2              Site_2;
  typedef typename K::Point_2             Point_2;

  typedef typename K::Comparison_result   result_type;

private:
  typedef typename K::Compare_x_2         Kernel_compare_x_2;

public:

  inline
  result_type operator()(const Point_2& p, const Point_2& q) const
  {
    return Kernel_compare_x_2()( p, q );
  }

  inline
  result_type operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    return Kernel_compare_x_2()( p.point(), q.point() );
  }
};

} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif //  CGAL_SEGMENT_DELAUNAY_GRAPH_2_COMPARE_X_2_H
