// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$

// $Id$
//
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_PREDICATES_ON_WEIGHTED_POINTS_CARTESIAN_3
#define CGAL_PREDICATES_ON_WEIGHTED_POINTS_CARTESIAN_3

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

namespace CGAL {

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
   return CGAL_NTS compare(dqx*dqx + dqy*dqy + dqz*dqz - qw,
                           drx*drx + dry*dry + drz*drz - rw);
}


//return the sign of the power test of weighted point (sx,sy,sz,sw)
//with respect to the smallest sphere orthogonal to
//p,q,r
template< class FT >
CGAL_MEDIUM_INLINE
Sign
in_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &sx, const FT &sy, const FT &sz, const FT  &sw)
{
 // Translate p to origin and compute determinants
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

  FT den = determinant(qq,qr,
			     qr,rr);
  FT detq = determinant(qpw,qr,
			      rpw,rr);
  FT detr = determinant(qq,qpw,
			      qr,rpw);

  // Smallest  smallest orthogonal sphere center
  // c =  detq/2*den  q + detr/2*den  r (origin at p)
  // square radius  c^2 - pw

  FT spx = sx-px;
  FT spy = sy-py;
  FT spz = sz-pz;
  FT ss = CGAL_NTS square(spx) + CGAL_NTS square(spy) + CGAL_NTS  square(spz);
  FT sq = spx*qpx + spy*qpy + spz*qpz;
  FT sr = spx*rpx + spy*rpy + spz*rpz;

  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );
  // return  sign of (c- s)^2 - (c^2 - pw) - sw    note that den >= 0 -
  return CGAL_NTS sign( den*(ss - sw + pw)- detq*sq - detr*sr);
}



// return the sign of the power test of weighted point (rx,ry,rz,rw)
 // with respect to the smallest sphere orthogoanal to
// p,q
template< class FT >
Sign
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
  return CGAL_NTS sign (dr2 - dp2/FT4 + dpr*dpw/dp2 - drw );
}




// return ON_UNBOUNDED_SIDE, ON_BOUNDARY or ON_BOUNDED_SIDE according
// to the position of the weighted circumcenter of pqrs
// with respect with the tertraedron formed by bare points in p, q,r,s
template <class FT>
Bounded_side
does_simplex_intersect_weighted_dual_supportC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw)
{
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

  FT den = determinant(qq,qr,qs,
			     qr,rr,rs,
			     qs,rs,ss);
  FT detq = determinant(qpw,qr,qs,
			      rpw,rr,rs,
			      spw,rs,ss);
  FT detr = determinant(qq,qpw,qs,
			      qr,rpw,rs,
			      qs,spw,ss);
  FT dets = determinant(qq,qr,qpw,
			      qr,rr,rpw,
			      qs,rs,spw);
  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );


  // The barycentrique coordinate of the smallest orthogonal sphere center
  // are  detq/2*den detr/2*den dets/2*den
  // and  1-(detq+ detr+dets)/2*den

  CGAL::Sign  sign1 = CGAL_NTS sign(FT(2)*den - detq -detr -dets);
  if (
   (CGAL_NTS sign(detq) == CGAL_NTS sign(den) || CGAL_NTS sign(detq)== ZERO) &&
   (CGAL_NTS sign(detr) == CGAL_NTS sign(den) || CGAL_NTS sign(detr)== ZERO) &&
   (CGAL_NTS sign(dets) == CGAL_NTS sign(den) || CGAL_NTS sign(dets)== ZERO) &&
   ( sign1 == POSITIVE || sign1 == ZERO )) { // inside or on boundary
    if (CGAL_NTS sign(detq) != ZERO &&
	CGAL_NTS sign(detr) != ZERO &&
	CGAL_NTS sign(dets) != ZERO &&
	sign1 != ZERO)
      return ON_BOUNDED_SIDE;
    else return ON_BOUNDARY ;
  }
  return ON_UNBOUNDED_SIDE;
}

// return ON_UNBOUNDED_SIDE, ON_BOUNDARY or ON_BOUNDED_SIDE according
// to the position of the  intersection if the radialaxis
// of weighted circumcenter of pqr with affine hull of bare p,q,r
// with respect with the triangle  formed by bare points in p, q,r.
template <class FT>
Bounded_side
does_simplex_intersect_weighted_dual_supportC3(
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

  FT den = determinant(qq,qr,
			     qr,rr);
  FT detq = determinant(qpw,qr,
			      rpw,rr);
  FT detr = determinant(qq,qpw,
			      qr,rpw);

  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );

   // The barycentrique coordinate of the smallest orthogonal sphere center
  // are  detq/2*den detr/2*den
  // and  1-(detq+ detr)/2*den

 CGAL::Sign  sign1 = CGAL_NTS sign(FT(2)*den - detq - detr);
 if (
   (CGAL_NTS sign(detq) == CGAL_NTS sign(den) || CGAL_NTS sign(detq)== ZERO) &&
   (CGAL_NTS sign(detr) == CGAL_NTS sign(den) || CGAL_NTS sign(detr)== ZERO) &&
   ( sign1 == POSITIVE || sign1 == ZERO )) { // inside or on boundary
   if ( CGAL_NTS sign(detq) != ZERO &&
	CGAL_NTS sign(detr) != ZERO &&
	sign1 != ZERO)
     return ON_BOUNDED_SIDE;
   else return ON_BOUNDARY;
 }
 return ON_UNBOUNDED_SIDE;
}


template <class FT>
Bounded_side
does_simplex_intersect_weighted_dual_supportC3(
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

  CGAL::Sign  sign1 = CGAL_NTS sign( qq + dw);
  CGAL::Sign  sign2 = CGAL_NTS sign( qq - dw);

  if( ( sign1 == POSITIVE || sign1 == ZERO ) &&
      ( sign2 == POSITIVE || sign2 == ZERO )) { // inside or on boundary
    if (sign1 != ZERO && sign2 != ZERO) return ON_BOUNDED_SIDE;
    else return ON_BOUNDARY;
  }
  return  ON_UNBOUNDED_SIDE;
}

//-------------------------------------------------------------------
} //namespace CGAL
//-------------------------------------------------------------------

#endif //CGAL_PREDICATES_ON_WEIGHTED_POINTS_CARTESIAN_3
