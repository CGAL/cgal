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
// file          : include/CGAL/smallest_radius_2.h
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_SMALLEST_RADIUS_2_H
#define CGAL_SMALLEST_RADIUS_2_H

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/smallest_radiusH2.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/smallest_radiusC2.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class R >
inline
R_FT_return(R)
squared_radius_smallest_circumcircle(const Point_2<R> &p,
				     const Point_2<R> &q) 
{
  
  Vector_2<R> v(p - q);
  return R_FT_return(R)((v*v)/4);
}

template <class R >
inline
R_FT_return(R)
squared_radius_smallest_circumcircle(const Point_2<R> &p,
				     const Point_2<R> &q,
				     const Point_2<R> &r) 
{
  typedef typename R::Point_2_base Point_two; 
  // compute the smallest radius directly
  return squared_radius_smallest_circumcircle((const Point_two&) p,
							      (const Point_two&) q,
							      (const Point_two&) r);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif // SMALLEST_RADIUS_2_H
