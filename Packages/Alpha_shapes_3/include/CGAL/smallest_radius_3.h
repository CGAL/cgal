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
// file          : include/CGAL/smallest_radius_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef SMALLEST_RADIUS_3_H
#define SMALLEST_RADIUS_3_H

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/smallest_radiusH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/smallest_radiusC3.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Point_3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class R >
typename R::FT
inline squared_radius_smallest_circumsphere(const Point_3<R> &p,
					    const Point_3<R> &q)
{
  Vector_3<R> v(p - q);
  return typename R::FT ((v*v)/R::FT(4));
}

template <class R >
typename R::FT
inline squared_radius_smallest_circumsphere(const Point_3<R> &p,
					    const Point_3<R> &q,
					    const Point_3<R> &r)
{
  typedef typename R::Point_3_base Point_three;
  return squared_radius_smallest_circumsphere((const Point_three&)p,
					      (const Point_three&)q,
					      (const Point_three&)r);
}

template <class R >
typename R::FT
inline squared_radius_circumsphere(const Point_3<R> &p,
				   const Point_3<R> &q,
				   const Point_3<R> &r,
				   const Point_3<R> &s)
{
  typedef typename R::Point_3_base Point_three;
  return squared_radius_circumsphere((const Point_three&)p,
				     (const Point_three&)q,
				     (const Point_three&)r,
				     (const Point_three&)s);  
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif // SMALLEST_RADIUS_3_H
