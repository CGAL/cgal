// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/side_of_smallest_circle_2.h
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// Additional predicates for predicates_on_points_2.h needed for computing
// alpha-shapes.
// ======================================================================

#ifndef CGAL_SIDE_OF_SMALLEST_CIRCLE_2_PLUS_H
#define CGAL_SIDE_OF_SMALLEST_CIRCLE_2_PLUS_H
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/predicates_on_points_2.h>

#include <CGAL/smallest_radius_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class R >
inline Bounded_side side_of_bounded_circle(const Point_2<R>& p,
					   const Point_2<R>& q,
					   const Point_2<R>& test) 
{
  
  Vector_2<R> v(test - (p + (q - p)/2));
  typename R::FT squared_distance = v*v;

  typename R::FT squared_radius = squared_radius_smallest_circumcircle(p, q);
 
  return (squared_radius > squared_distance) ?
    ON_BOUNDED_SIDE :
    ((squared_radius < squared_distance) ?
     ON_UNBOUNDED_SIDE :
     ON_BOUNDARY);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif  // SIDE_OF_SMALLEST_CIRCLE_2_PLUS_H
