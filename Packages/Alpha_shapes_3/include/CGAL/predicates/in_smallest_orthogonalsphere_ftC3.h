// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$

// $Revision$ $Date$
// $Name$
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

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


template< class FT >
Bounded_side
in_smallest_orthogonalsphereC3(
 const FT &px, const FT &py, const FT &pz, const FT  &pw,
 const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
 const FT &rx, const FT &ry, const FT &rz, const FT  &rw)
{
  FT FT2(2);
  FT FT4(4);
  FT dpx = px - qx;
  FT dpy = py - qy;
  FT dpz = pz - qz;
  FT dpw = pw - qw;
  FT dp2 = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) + CGAL_NTS square(dpz); 
  FT drx = rx - (px + qx)/FT2;
  FT dry = ry - (py + qy)/FT2;
  FT drz = rz - (pz + qz)/FT2;
  FT drw = rw - (pw + qw)/FT2;
  FT dr2 = CGAL_NTS square(drx) + CGAL_NTS square(dry) + CGAL_NTS square(drz); 
  FT dpr = dpx*drx + dpy*dry +dpz*drz;
  return Bounded_side(CGAL_NTS sign (dr2 - dp2/FT4 + dpr*dpw/dp2 + drw ));
}
//------------------------------------------------------------------- 
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONALSPHEREC3_H
