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
// file          : include/CGAL/in_smallest_orthogonalcircleC2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IN_SMALLEST_ORTHOGONALCIRCLEC2_H 
#define CGAL_IN_SMALLEST_ORTHOGONALCIRCLEC2_H

#ifndef CGAL_POINTC2_H
#include <CGAL/Cartesian/Point_2.h>
#endif // CGAL_POINTC2_H

#ifndef CGAL_WEIGHTED_POINT_H
#include <CGAL/Weighted_point.h>
#endif // CGAL_WEIGHTED_POINT_H

#include <CGAL/predicates/in_smallest_orthogonalcircle_ftC2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
in_smallest_orthogonalcircle(
		   const Weighted_point<Point_2< Cartesian<FT> >, FT > &p,
		   const Weighted_point<Point_2< Cartesian<FT> >, FT > &q,
		   const Weighted_point<Point_2< Cartesian<FT> >, FT > &t) 
{

 
  FT px(p.point().x());
  FT py(p.point().y());
  FT pw(p.weight());
  FT qx(q.point().x());
  FT qy(q.point().y());
  FT qw(q.weight());
  FT tx(t.point().x());
  FT ty(t.point().y());
  FT tw(t.weight());

  return in_smallest_orthogonalcircleC2(px, py, pw,
					qx, qy, qw,
					tx, ty, tw);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONALCIRCLEC2_H
