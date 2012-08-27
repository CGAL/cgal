// Copyright (c) 2004   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H
#define CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H

namespace CGAL {

template <class FT>
void
determinants_for_weighted_circumcenterC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                FT &num_x,  FT &num_y, FT &num_z, FT& den)
{
  // translate origin to p
  // and compute determinants for weighted_circumcenter and
  // circumradius
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz) - qw + pw;
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) +
           CGAL_NTS square(rpz) - rw + pw;
  FT spx = sx-px;
  FT spy = sy-py;
  FT spz = sz-pz;
  FT sp2 = CGAL_NTS square(spx) + CGAL_NTS square(spy) +
           CGAL_NTS square(spz) - sw + pw;

  num_x = determinant(qpy,qpz,qp2,
			    rpy,rpz,rp2,
			    spy,spz,sp2);
  num_y = determinant(qpx,qpz,qp2,
			    rpx,rpz,rp2,
			    spx,spz,sp2);
  num_z = determinant(qpx,qpy,qp2,
			    rpx,rpy,rp2,
			    spx,spy,sp2);
  den   = determinant(qpx,qpy,qpz,
			    rpx,rpy,rpz,
			    spx,spy,spz);
}


template < class FT>
void
weighted_circumcenterC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                FT &x, FT &y, FT &z)
{
  // this function  compute the weighted circumcenter point only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   num_x,  num_y, num_z,den);

  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );
   FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
}

template < class FT>
void
weighted_circumcenterC3(
		const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                FT &x, FT &y, FT &z, FT &w)
{
  // this function  compute the weighted circumcenter point
  // and the squared weighted circumradius

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   num_x,  num_y, num_z, den);

  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );
   FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;

  w = (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
      *CGAL_NTS square(inv) - pw;
}


template< class FT >
FT
squared_radius_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw,
  const FT &sx, const FT &sy, const FT &sz, const FT  &sw)
{

  // this function  compute the squared weighted circumradius only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   num_x, num_y, num_z,den);

  CGAL_triangulation_assertion( ! CGAL_NTS is_zero(den) );
   FT inv = FT(1)/(FT(2) * den);

   return
    (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
    *CGAL_NTS square(inv) - pw;
}


template <class FT>
void
determinants_for_weighted_circumcenterC3(
	        const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		FT &num_x,  FT &num_y, FT &num_z, FT &den)
{
  // translate origin to p
  // and compute determinants for weighted_circumcenter and
  // circumradius

  // Translate s to origin to simplify the expression.
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz) - qw + pw;
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) +
           CGAL_NTS square(rpz) - rw + pw;

  FT sx = qpy*rpz-qpz*rpy;
  FT sy = qpz*rpx-qpx*rpz;
  FT sz = qpx*rpy-qpy*rpx;

  // The following determinants can be developped and simplified.
  //
  // FT num_x = determinant(qpy,qpz,qp2,
  //                              rpy,rpz,rp2,
  //                              sy,sz,FT(0));
  // FT num_y = determinant(qpx,qpz,qp2,
  //                              rpx,rpz,rp2,
  //                              sx,sz,FT(0));
  // FT num_z = determinant(qpx,qpy,qp2,
  //                              rpx,rpy,rp2,
  //                              sx,sy,FT(0));

  num_x = qp2 * determinant(rpy,rpz,sy,sz)
        - rp2 * determinant(qpy,qpz,sy,sz);

  num_y = qp2 * determinant(rpx,rpz,sx,sz)
	- rp2 * determinant(qpx,qpz,sx,sz);

  num_z = qp2 * determinant(rpx,rpy,sx,sy)
	- rp2 * determinant(qpx,qpy,sx,sy);

  den   = determinant(qpx,qpy,qpz,
			    rpx,rpy,rpz,
			    sx,sy,sz);
}

template < class FT >
void
weighted_circumcenterC3(
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		  FT &x, FT &y, FT &z)
{
  // this function  compute the weighted circumcenter point only

// Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   num_x,  num_y, num_z, den);

  CGAL_triangulation_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
}


template < class FT >
void
weighted_circumcenterC3(
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  const FT &rx, const FT &ry, const FT &rz, const FT &rw,
		  FT &x, FT &y, FT &z, FT &w)
{
  // this function  compute the weighted circumcenter and
  // the weighted squared circumradius

// Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   num_x,  num_y, num_z, den);

  CGAL_triangulation_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;

  w = (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
      *CGAL_NTS square(inv)  - pw;
}


template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw,
  const FT &rx, const FT &ry, const FT &rz, const FT  &rw)
{
  // this function  compute the weighted squared circumradius only

// Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   num_x,  num_y, num_z, den);

  CGAL_triangulation_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

 return
    (CGAL_NTS square(num_x)+CGAL_NTS square(num_y)+CGAL_NTS square(num_z))
     *CGAL_NTS square(inv)  - pw;
}



template < class FT >
void
weighted_circumcenterC3(
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  FT &x, FT &y, FT &z)
{
// this function  compute the weighted circumcenter point only
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT inv = FT(1)/(FT(2)*qp2);
  FT alpha = 1/FT(2) + (pw-qw)*inv;

  x = px + alpha * qpx;
  y = py + alpha * qpy;
  z = pz + alpha * qpz;
}


template < class FT >
void
weighted_circumcenterC3(
                  const FT &px, const FT &py, const FT &pz, const FT &pw,
		  const FT &qx, const FT &qy, const FT &qz, const FT &qw,
		  FT &x, FT &y, FT &z, FT &w)
{
 // this function  compute the weighted circumcenter point and
  // the weighted circumradius
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT inv = FT(1)/(FT(2)*qp2);
  FT alpha = 1/FT(2) + (pw-qw)*inv;

  x = px + alpha * qpx;
  y = py + alpha * qpy;
  z = pz + alpha * qpz;

  w = CGAL_NTS square(alpha)*qp2 - pw;
}


template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{
  // this function  computes
  // the weighted circumradius only
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT inv = FT(1)/(FT(2)*qp2);
  FT alpha = 1/FT(2) + (pw-qw)*inv;

  return  CGAL_NTS square(alpha)*qp2 - pw;
}


template< class FT >
FT
power_productC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{
  // computes the power product of two weighted points
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  return qp2 - pw - qw ;
}

template < class RT , class We>
void
radical_axisC3(const RT &px, const RT &py, const RT &pz, const We & /* pw */,
	       const RT &qx, const RT &qy, const RT &qz, const We & /* qw */,
	       const RT &rx, const RT &ry, const RT &rz, const We & /* rw */,
	       RT &a, RT &b, RT& c )
{
  RT dqx=qx-px, dqy=qy-py, dqz=qz-pz, drx=rx-px, dry=ry-py, drz=rz-pz;

  //il manque des tests...

  a= RT(1)*determinant(dqy, dqz, dry, drz);
  b= - RT(1)*determinant(dqx, dqz, drx, drz);
  c= RT(1)*determinant(dqx, dqy, drx, dry);
}

// function used in critical_squared_radiusC3
// power ( t, tw) with respect to
// circle orthogonal (p,pw), (q,qw), (r,rw), (s,sw)
template < class FT>
FT
power_to_orthogonal_sphereC3(
                const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                const FT &tx, const FT &ty, const FT &tz, const FT &tw)
{
   //to get the value of the determinant
   // We translate the points so that t  becomes the origin.
    FT dpx = px - tx;
    FT dpy = py - ty;
    FT dpz = pz - tz;
    FT dpt = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) +
             CGAL_NTS square(dpz) - pw + tw ;
    FT dqx = qx - tx;
    FT dqy = qy - ty;
    FT dqz = qz - tz;
    FT dqt = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) +
             CGAL_NTS square(dqz) - qw + tw;
    FT drx = rx - tx;
    FT dry = ry - ty;
    FT drz = rz - tz;
    FT drt = CGAL_NTS square(drx) + CGAL_NTS square(dry) +
             CGAL_NTS square(drz) - rw + tw;
    FT dsx = sx - tx;
    FT dsy = sy - ty;
    FT dsz = sz - tz;
    FT dst = CGAL_NTS square(dsx) + CGAL_NTS square(dsy) +
             CGAL_NTS square(dsz) - sw + tw;

    return determinant(dpx, dpy, dpz, dpt,
			     dqx, dqy, dqz, dqt,
			     drx, dry, drz, drt,
			     dsx, dsy, dsz, dst);

}



// compute the critical weight tw
// where weighted point t is orthogonal to weighted points p, q,r,s
template < class FT>
FT
critical_squared_radiusC3(
                const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                const FT &tx, const FT &ty, const FT &tz, const FT &  )
{
  // the 5x5 det  is a linear function of tw ff(tw)= ff(0) + tw ff(1)
  // the critical value for tw is  - ff(0)/( ff(1) - ff(0))


    FT ff0 =  power_to_orthogonal_sphereC3(px, py, pz, pw,
					   qx, qy, qz, qw,
					   rx, ry, rz, rw,
					   sx, sy, sz, sw,
					   tx, ty, tz, FT(0));

    FT ff1 = power_to_orthogonal_sphereC3(px, py, pz, pw,
					    qx, qy, qz, qw,
					    rx, ry, rz, rw,
					    sx, sy, sz, sw,
					    tx, ty, tz, FT(1));

    return   -ff0/(ff1 - ff0);
}



 // I will use this to test if the radial axis of three spheres
  // intersect the triangle formed by the centers.
//   // resolution of the system (where we note c the center)
//   // |       dc^2 = cw + rw
//   // |  (dp-dc)^2 = pw + cw
//   // |  (dq-dc)^2 = qw + cw
//   // |         dc = Lamdba*dp + Mu*dq

//   FT FT2(2);
//   FT dpx = px-rx;
//   FT dpy = py-ry;
//   FT dpz = pz-rz;
//   FT dp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)+CGAL_NTS  square(dpz);
//   FT dpp = dp-pw+rw;
//   FT dqx = qx-rx;
//   FT dqy = qy-ry;
//   FT dqz = qz-rz;
//   FT dq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)+CGAL_NTS square(dqz);
//   FT dqq = dq-qw+rw;
//   FT dpdq = dpx*dqx+dpy*dqy+dpz*dqz;
//   FT denom = FT2*(dp*dq-CGAL_NTS square(dpdq));
//   FT Lambda = (dpp*dq-dqq*dpdq)/denom;
//   FT Mu = (dqq*dp-dpp*dpdq)/denom;

//   return (CGAL_NTS square(Lambda)*dp+CGAL_NTS square(Mu)*dq
// 	  +FT2*Lambda*Mu*dpdq - rw);




} //namespace CGAL
#endif //CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_3_H
