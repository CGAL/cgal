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

#ifndef CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALSPHERE_FTC3_H 
#define CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALSPHERE_FTC3_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_orthogonalsphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &sx, const FT &sy, const FT &sz, const FT  &sw)
{
  FT FT4(4);
  FT dpx = px-sx;
  FT dpy = py-sy;
  FT dpz = pz-sz;
  FT dpp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)+CGAL_NTS square(dpz)
           -pw+sw;
  FT dqx = qx-sx;
  FT dqy = qy-sy;
  FT dqz = qz-sz;
  FT dqq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)+CGAL_NTS square(dqz)
           -qw+sw;
  FT drx = rx-sx;
  FT dry = ry-sy;
  FT drz = rz-sz;
  FT drr = CGAL_NTS square(drx)+CGAL_NTS square(dry)+CGAL_NTS square(drz)
           -rw+sw;

  FT det0 = det3x3_by_formula(dpx,dpy,dpz,dqx,dqy,dqz,drx,dry,drz);
  
  FT det1 = det3x3_by_formula(dpp,dpy,dpz,dqq,dqy,dqz,drr,dry,drz);
  FT det2 = det3x3_by_formula(dpx,dpp,dpz,dqx,dqq,dqz,drx,drr,drz);
  FT det3 = det3x3_by_formula(dpx,dpy,dpp,dqx,dqy,dqq,drx,dry,drr);

  return
    (CGAL_NTS square(det1)+CGAL_NTS square(det2)+CGAL_NTS square(det3))/
    (FT4*CGAL_NTS square(det0)) - sw;
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonalsphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw)
{
  // resolution of the system (where we note c the center)
  // |       dc^2 = cw + rw
  // |  (dp-dc)^2 = pw + cw
  // |  (dq-dc)^2 = qw + cw
  // |         dc = Lamdba*dp + Mu*dq

  FT FT2(2);
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
  FT denom = FT2*(dp*dq-CGAL_NTS square(dpdq));
  FT Lambda = (dpp*dq-dqq*dpdq)/denom;
  FT Mu = (dqq*dp-dpp*dpdq)/denom;

  return (CGAL_NTS square(Lambda)*dp+CGAL_NTS square(Mu)*dq
	  +FT2*Lambda*Mu*dpdq - rw);
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonalsphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{ 
  FT FT4(4);
  FT dp = CGAL_NTS square(px-qx)+CGAL_NTS square(py-qy)+CGAL_NTS square(pz-qz);

  return (CGAL_NTS square(dp-pw+qw)/(FT4*dp)-qw);
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALSPHERE_ftC3_H
