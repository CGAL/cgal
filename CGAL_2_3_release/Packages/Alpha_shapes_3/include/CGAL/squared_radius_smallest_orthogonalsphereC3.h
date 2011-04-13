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
// file          : include/CGAL/squared_radius_smallest_orthogonalsphereC3.h
// package       : Alpha_shapes_3 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALSPHEREC3_H 
#define CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALSPHEREC3_H

#ifndef CGAL_POINTC3_H
#include <CGAL/Cartesian/Point_3.h>
#endif // CGAL_POINTC3_H

#ifndef CGAL_WEIGHTED_POINT_H
#include <CGAL/Weighted_point.h>
#endif // CGAL_WEIGHTED_POINT_H

#include <CGAL/constructions/squared_radius_smallest_orthogonalsphere_ftC3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT 
squared_radius_orthogonalsphere(
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &p,
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &q,
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &r,
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &s) 
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
  FT sx(s.point().x());
  FT sy(s.point().y());
  FT sz(s.point().z());
  FT sw(s.weight());
  return squared_radius_orthogonalsphereC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw);
}

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT 
squared_radius_smallest_orthogonalsphere(
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &p,
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &q,
  const Weighted_point<Point_3< Cartesian<FT> >, FT> &r) 
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

  return squared_radius_smallest_orthogonalsphereC3(px, py, pz, pw,
						    qx, qy, qz, qw,
						    rx, ry, rz, rw);
}

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT 
squared_radius_smallest_orthogonalsphere(
   const Weighted_point<Point_3< Cartesian<FT> >, FT> &p,
   const Weighted_point<Point_3< Cartesian<FT> >, FT> &q) 
{
  FT px(p.point().x());
  FT py(p.point().y());
  FT pz(p.point().z());
  FT pw(p.weight());
  FT qx(q.point().x());
  FT qy(q.point().y());
  FT qz(q.point().z());
  FT qw(q.weight());

  return squared_radius_smallest_orthogonalsphereC3(px, py, pz, pw,
						    qx, qy, qz, qw);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALSPHERE_C3_H
