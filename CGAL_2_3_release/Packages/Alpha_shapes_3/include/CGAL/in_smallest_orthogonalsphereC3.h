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
// file          : include/CGAL/in_smallest_orthogonalsphereC3.h
// package       : Alpha_shapes_3 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IN_SMALLEST_ORTHOGONALSPHEREC3_H 
#define CGAL_IN_SMALLEST_ORTHOGONALSPHEREC3_H

#ifndef CGAL_POINTC3_H
#include <CGAL/Cartesian/Point_3.h>
#endif // CGAL_POINTC3_H

#ifndef CGAL_WEIGHTED_POINT_H
#include <CGAL/Weighted_point.h>
#endif // CGAL_WEIGHTED_POINT_H

#include <CGAL/predicates/in_smallest_orthogonalsphere_ftC3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
in_smallest_orthogonalsphere(
			     const Weighted_point<Point_3< Cartesian<FT> >, FT > &p,
			     const Weighted_point<Point_3< Cartesian<FT> >, FT > &q, 
			     const Weighted_point<Point_3< Cartesian<FT> >, FT > &r,
			     const Weighted_point<Point_3< Cartesian<FT> >, FT > &t) 
{

 
  FT px(p.point().x());
  FT py(p.point().y());
  FT pz(p.point().z());
  FT pw(p.weight());
  FT qx(q.point().x());
  FT qy(q.point().y());
  FT qz(q.point().z());
  FT qw(q.weight());
  FT rx(r.point().x());
  FT ry(r.point().y());
  FT rz(r.point().z());
  FT rw(r.weight());
  FT tx(t.point().x());
  FT ty(t.point().y());
  FT tz(t.point().z());
  FT tw(t.weight());

  return in_smallest_orthogonalsphereC3(px, py, pz, pw,
					qx, qy, qz, qw,
					rx, ry, rz, rw,
					tx, ty, tz, tw);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONALSPHEREC3_H
