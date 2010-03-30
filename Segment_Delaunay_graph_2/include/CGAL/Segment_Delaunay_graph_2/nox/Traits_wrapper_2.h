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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_TRAITS_WRAPPER_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_TRAITS_WRAPPER_2_H

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

CGAL_BEGIN_NAMESPACE

template<class Gt_base>
class Segment_Delaunay_graph_traits_wrapper_2 : public Gt_base
{
private:
  typedef Segment_Delaunay_graph_traits_wrapper_2<Gt_base> Self;

public:
  struct Triangle_2 {};

  Segment_Delaunay_graph_traits_wrapper_2() {}

  Segment_Delaunay_graph_traits_wrapper_2(const Gt_base&) {}

  struct Side_of_oriented_circle_2
  {
    typedef typename Gt_base::Point_2 Point_2;

    inline
    Sign operator()(const Point_2& p1, const Point_2& p2,
		    const Point_2& p3, const Point_2& q) const
    {
      return CGAL::side_of_oriented_circle(p1, p2, p3, q);
    }
  };

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2();
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_NOX_TRAITS_WRAPPER_2_H
