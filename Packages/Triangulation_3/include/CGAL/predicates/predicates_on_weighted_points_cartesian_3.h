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

#ifndef CGAL_IN_SMALLEST_ORTHOGONAL_SPHERE_FTC3_H 
#define CGAL_IN_SMALLEST_ORTHOGONAL_SPHERE_FTC3_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class FT >
Comparison_result
compare_power_distanceC3(
		  const FT &px, const FT &py, const FT &pz, 
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  const FT &rx, const FT &ry, const FT &rz, const FT &rw)
{
 FT dqx = qx - px;
 FT dqy = qy - py;
 FT dqz = qz - pz;
 FT drx = rx - px;
 FT dry = ry - py;
 FT drz = rz - pz;
   return Comparison_result(CGAL_NTS sign (
      (dqx*dqx + dqy*dqy + dqz*dqz - qw )
    - (drx*drx + dry*dry + drz*drz - rw ) ));
}

template< class FT >
CGAL_MEDIUM_INLINE
Bounded_side
in_smallest_orthogonal_sphereC3(
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
in_smallest_orthogonal_sphereC3(
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

template <class FT>
bool
does_affine_dual_intersectC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw)
{
  // returns true if weighted circumcenter of pqrs is in
  // tetrahedron pqrs or on boundary
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  
  FT spx = sx-px;
  FT spy = sy-py;
  FT spz = sz-pz;

  FT qq = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
  FT rr = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + CGAL_NTS square(rpz);
  FT ss = CGAL_NTS square(spx) + CGAL_NTS square(spy) + CGAL_NTS square(spz);
  FT qr = qpx*rpx + qpy*rpy + qpz*rpz;
  FT rs = rpx*spx + rpy*spy + rpz*spz; 
  FT qs = spx*qpx + spy*qpy + spz*qpz;
  FT qpw = qq - qw + pw ;
  FT rpw = rr - rw + pw ;
  FT spw = ss - sw + pw ;
 
  FT den = det3x3_by_formula(qq,qr,qs,
			     qr,rr,rs,
			     qs,rs,ss);
  FT detq = det3x3_by_formula(qpw,qr,qs,
			      rpw,rr,rs,
			      spw,rs,ss);
  FT detr = det3x3_by_formula(qq,qpw,qs,
			      qr,rpw,rs,
			      qs,spw,ss);
  FT dets = det3x3_by_formula(qq,qr,qpw,
			      qr,rr,rpw,
			      qs,rs,spw);
  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );

 CGAL::Sign  sign = sign(FT(2)*den - detq -detr -dets)
  return  
   (CGAL_NTS sign(detq) == CGAL_NTS sign(den) || CGAL_NTS sign(detq)== ZERO) &&
   (CGAL_NTS sign(detr) == CGAL_NTS sign(den) || CGAL_NTS sign(detr)== ZERO) &&
   (CGAL_NTS sign(dets) == CGAL_NTS sign(den) || CGAL_NTS sign(dets)== ZERO) &&
   ( sign == POSITIVE || sign == ZERO );
}

template <class FT>
bool
does_affine_dual_intersectC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw)
{
  // returns true if weighted circumcenter of pqr is in
  // triangle pqr or on boundary
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  
  FT qq = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
  FT rr = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + CGAL_NTS square(rpz);
  FT qr = qpx*rpx + qpy*rpy + qpz*rpz;

  FT qpw = qq - qw + pw ;
  FT rpw = rr - rw + pw ;
  
  FT den = det2,2_by_formula(qq,qr,
			     qr,rr);
  FT detq = det3x3_by_formula(qpw,qr,
			      rpw,rr);
  FT detr = det3x3_by_formula(qq,qpw,
			      qr,rpw);

  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );

 CGAL::Sign  sign = CGAL_NTS sign(FT(2)*den - detq - detr)
  return  
   (CGAL_NTS sign(detq) == CGAL_NTS sign(den) || CGAL_NTS sign(detq)== ZERO) &&
   (CGAL_NTS sign(detr) == CGAL_NTS sign(den) || CGAL_NTS sign(detr)== ZERO) &&
   ( sign == POSITIVE || sign == ZERO );
}


template <class FT>
bool
does_affine_dual_intersectC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw)
{
  // returns true if weighted circumcenter of pq is in
  // segment  pq or on boundary
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  
  FT qq = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
  FT dw = pw - qw;
  
  CGAL::Sign  sign1 = CGAL_NTS sign( qq - dw);
  CGAL::Sign  sign2 = CGAL_NTS sign( qq - dw);
   
  return 
    ( sign1 == POSITIVE || sign == ZERO ) &&
    ( sign2 == POSITIVE || sign == ZERO );
}

//------------------------------------------------------------------- 
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONAL_SPHEREC3_H
