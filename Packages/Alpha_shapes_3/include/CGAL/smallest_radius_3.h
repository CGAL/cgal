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
// release       :
// release_date  :
//
// file          : include/CGAL/smallest_radius_3.h
// package       : Alpha_shapes_3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_SMALLEST_RADIUS_3_H
#define CGAL_SMALLEST_RADIUS_3_H

#include <CGAL/Point_3.h>

CGAL_BEGIN_NAMESPACE

template <class R >
CGAL_KERNEL_MEDIUM_INLINE
typename R::FT
squared_radius_smallest_circumsphere(const Point_3<R> &p,
				     const Point_3<R> &q)
{
  Vector_3<R> v(p - q);
  return typename R::FT ((v*v)/R::FT(4));
}

CGAL_END_NAMESPACE

#endif // CGAL_SMALLEST_RADIUS_3_H
