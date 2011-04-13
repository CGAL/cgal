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
// file          : include/CGAL/predicates/in_smallest_orthogonalsphere_ftC3.h
// package       : Alpha_shapes_3 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IN_SMALLEST_ORTHOGONALSPHERE_FTC3_H 
#define CGAL_IN_SMALLEST_ORTHOGONALSPHERE_FTC3_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template< class FT >
CGAL_MEDIUM_INLINE
Bounded_side
in_smallest_orthogonalsphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &tx, const FT &ty, const FT &tz, const FT  &tw)
{
  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dpz = pz-rz;
  FT dp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)+CGAL_NTS square(dpz);
  FT dpp = dp-pw+rw;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dqz = qz-rz;
  FT dq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)+CGAL_NTS square(dqz);
  FT dqq = dq-qw+rw;
  FT dpdq = dpx*dqx+dpy*dqy+dpz*dqz;
  FT denom = dp*dq-CGAL_NTS square(dpdq);
  FT Lambda = (dpp*dq-dqq*dpdq)/denom;
  FT Mu = (dqq*dp-dpp*dpdq)/denom;
  
  FT dtx = tx-rx;
  FT dty = ty-ry;
  FT dtz = tz-rz;
  
  return Bounded_side 
    (CGAL_NTS sign(
	-(CGAL_NTS square(dtx)+CGAL_NTS square(dty)+CGAL_NTS square(dtz)-tw+rw)
        +Lambda*(dpx*dtx+dpy*dty+dpz*dtz)+Mu*(dqx*dtx+dqy*dty+dqz*dtz)));
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONALSPHEREC3_H
