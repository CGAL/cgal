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
// file          : include/CGAL/side_of_smallest_sphere_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

/* Additional predicates for predicates_on_points_3.h needed for computing
   alpha-shapes.                                                          */

#ifndef SIDE_OF_SMALLEST_SPHERE_3_H
#define SIDE_OF_SMALLEST_SPHERE_3_H
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/predicates_on_points_3.h>
#include <CGAL/smallest_radius_3.h>

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/side_of_smallest_sphereH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/side_of_smallest_sphereC3.h>
#endif // CGAL_CARTESIAN_H

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class R >
inline Bounded_side side_of_bounded_sphere(const Point_3<R> &p,
					   const Point_3<R> &q,                                              
					   const Point_3<R> &test) 
{
  Vector_3<R> v(test - (p + (q - p)/2));
  typename R::FT squared_distance = v*v;

  typename R::FT squared_radius = squared_radius_smallest_circumsphere(p, q);
 
  return (squared_radius > squared_distance) ?
                      ON_BOUNDED_SIDE :
                      ((squared_radius < squared_distance) ?
                                            ON_UNBOUNDED_SIDE :
                                            ON_BOUNDARY);
}

template <class R >
inline Bounded_side side_of_bounded_sphere(const Point_3<R> &p,
					   const Point_3<R> &q,  
					   const Point_3<R> &r,
					   const Point_3<R> &test) 
{
  typedef typename R::Point_3_base Point_three;
  return side_of_bounded_sphere((const Point_three&)p,
				(const Point_three&)q,
				(const Point_three&)r,
				(const Point_three&)test);
  
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif  // SIDE_OF_SMALLEST_SPHERE_3_H
