// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          :
//        include/CGAL/predicates/Segment_Voronoi_diagram_predicates_C2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_PREDICATES_2_H

#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_C2.h>
#include <CGAL/predicates/Svd_are_same_points_C2.h>
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
