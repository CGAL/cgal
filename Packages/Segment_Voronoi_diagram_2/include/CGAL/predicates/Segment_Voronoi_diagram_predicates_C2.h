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




#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H

#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_C2.h>
#include <CGAL/predicates/Svd_are_same_points_C2.h>
#include <CGAL/predicates/Svd_oriented_side_of_bisector_C2.h>
#include <CGAL/predicates/Svd_incircle_2.h>
#include <CGAL/predicates/Svd_finite_edge_interior_2.h>
#include <CGAL/predicates/Svd_infinite_edge_interior_2.h>
#include <CGAL/predicates/Svd_is_degenerate_edge_2.h>
#include <CGAL/predicates/Svd_do_intersect_C2.h>
#include <CGAL/predicates/Svd_oriented_side_C2.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#if 1
template< class K >
class Svd_are_same_points_2
{
public:
  typedef typename K::Point_2      Point_2;

  typedef typename K::Compare_x_2  compare_x_2;
  typedef typename K::Compare_y_2  compare_y_2;
public:

  bool operator()(const Point_2& p, const Point_2& q) const
  {
    if ( compare_x_2()(p, q) != EQUAL ) { return false; }
    Comparison_result res = compare_y_2()(p, q);
    return res == EQUAL;
  }
};
#endif

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template< class K >
class Svd_is_endpoint_of_segment_2
{
public:
  typedef typename K::Point_2                       Point_2;
  typedef typename K::Segment_2                     Segment_2;
  typedef typename CGAL::Svd_are_same_points_2<K>   Are_same_points_2;

  typedef typename K::Compare_x_2  compare_x_2;
  typedef typename K::Compare_y_2  compare_y_2;
public:

  bool operator()(const Point_2& p, const Segment_2& s) const
  {
    Are_same_points_2 are_same_points;
    return ( are_same_points(p, s.source()) ||
	     are_same_points(p, s.target()) );
  }
};


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H
