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
// file          : include/CGAL/constructions/squared_radius_smallest_orthogonalcircle_ftC2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALCIRCLE_FTC2_H 
#define CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALCIRCLE_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_orthogonalcircleC2(
  const FT &px, const FT &py, const FT  &pw,
  const FT &qx, const FT &qy, const FT  &qw,  
  const FT &rx, const FT &ry, const FT  &rw)

{
  FT FT4(4);
  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dpp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)-pw+rw;
  FT dqq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)-qw+rw;

  FT det0 = det2x2_by_formula(dpx, dpy, dqx, dqy);
  
  FT det1 = det2x2_by_formula(dpp, dpy, dqq, dqy);

  FT det2 = det2x2_by_formula(dpx, dpp, dqx, dqq);

  return 
    (CGAL_NTS square(det1)+CGAL_NTS square(det2))/
                                  (FT4*CGAL_NTS square(det0)) - rw;
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonalcircleC2(
  const FT &px, const FT &py, const FT  &pw,
  const FT &qx, const FT &qy, const FT  &qw)

{
  FT FT4(4);
  FT dpz = CGAL_NTS square(px-qx)+CGAL_NTS square(py-qy);

  return (CGAL_NTS square(dpz-pw+qw)/(FT4*dpz)-qw);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALCIRCLE_ftC2_H
