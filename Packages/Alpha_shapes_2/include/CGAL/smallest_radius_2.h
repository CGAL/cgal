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

#include <CGAL/Point_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class R >
CGAL_KERNEL_MEDIUM_INLINE
typename R::FT
squared_radius_smallest_circumcircle(const Point_2<R> &p,
				     const Point_2<R> &q) 
{
  typedef typename R::FT Return_type;
  Vector_2<R> v(p - q);
  return (Return_type ((v*v)/4));
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif // SMALLEST_RADIUS_2_H

