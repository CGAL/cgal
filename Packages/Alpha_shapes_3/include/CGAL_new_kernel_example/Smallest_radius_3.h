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
// file          : include/CGAL/Smallest_radius_3.h
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
#include <CGAL/Smallest_radiusH3.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Smallest_radiusC3.h>
#endif // CGAL_CARTESIAN_H

#include <CGAL/Point_3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class R >
R_FT_return(R)
inline squared_radius_smallest_circumsphere(const Point_3<R> &p,
					    const Point_3<R> &q)
{
  Vector_3<R> v(p - q);
  return R_FT_return(R)((v*v)/R::FT(4));
}

template <class R >
R_FT_return(R)
inline squared_radius_smallest_circumsphere(const Point_3<R> &p,
					    const Point_3<R> &q,
					    const Point_3<R> &r)
{
  return squared_radius_smallest_circumsphere((const R::Point_3&)p,
					      (const R::Point_3&)q,
					      (const R::Point_3&)r);
}

template <class R >
R_FT_return(R)
inline squared_radius_smallest_circumsphere(const Point_3<R> &p,
					    const Point_3<R> &q,
					    const Point_3<R> &r,
					    const Point_3<R> &s)
{
  return squared_radius_smallest_circumsphere((const R::Point_3&)p,
					      (const R::Point_3&)q,
					      (const R::Point_3&)r,
					      (const R::Point_3&)s);
  
}
//-------------------------------------------------------------------
//EXEMPLE D'ANDREAS
class Squared_radius_smallest_circumsphere
{
  public:
    template <class T>
    bool
    operator()(const T& p, const T& q, const T& r, const T& s) const
    { return squared_radius_smallest_circumsphere(p,q,r,s); }
};

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif // SMALLEST_RADIUS_3_H
