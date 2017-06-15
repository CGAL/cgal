// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Herve Bronnimann, Mariette Yvinec

#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
#define CGAL_CONSTRUCTIONS_KERNEL_FTC3_H

#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>

namespace CGAL {

template < class FT >
CGAL_KERNEL_INLINE
void
midpointC3( const FT &px, const FT &py, const FT &pz,
            const FT &qx, const FT &qy, const FT &qz,
            FT &x, FT &y, FT &z)
{
  x = (px + qx) / 2;
  y = (py + qy) / 2;
  z = (pz + qz) / 2;
}

template < class FT >
void
barycenterC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &p2z,
             FT &x, FT &y, FT &z)
{
   FT w2 = 1 - w1;
   x = w1 * p1x + w2 * p2x;
   y = w1 * p1y + w2 * p2y;
   z = w1 * p1z + w2 * p2z;
}

template < class FT >
void
barycenterC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &p2z, const FT &w2,
             FT &x, FT &y, FT &z)
{
   FT sum = w1 + w2;
   CGAL_kernel_assertion(sum != 0);
   x = (w1 * p1x + w2 * p2x) / sum;
   y = (w1 * p1y + w2 * p2y) / sum;
   z = (w1 * p1z + w2 * p2z) / sum;
}

template < class FT >
void
barycenterC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &p2z, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &p3z,
             FT &x, FT &y, FT &z)
{
   FT w3 = 1 - w1 - w2;
   x = w1 * p1x + w2 * p2x + w3 * p3x;
   y = w1 * p1y + w2 * p2y + w3 * p3y;
   z = w1 * p1z + w2 * p2z + w3 * p3z;
}

template < class FT >
void
barycenterC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &p2z, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &p3z, const FT &w3,
             FT &x, FT &y, FT &z)
{
   FT sum = w1 + w2 + w3;
   CGAL_kernel_assertion(sum != 0);
   x = (w1 * p1x + w2 * p2x + w3 * p3x) / sum;
   y = (w1 * p1y + w2 * p2y + w3 * p3y) / sum;
   z = (w1 * p1z + w2 * p2z + w3 * p3z) / sum;
}

template < class FT >
void
barycenterC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &p2z, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &p3z, const FT &w3,
             const FT &p4x, const FT &p4y, const FT &p4z,
             FT &x, FT &y, FT &z)
{
   FT w4 = 1 - w1 - w2 - w3;
   x = w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x;
   y = w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y;
   z = w1 * p1z + w2 * p2z + w3 * p3z + w4 * p4z;
}

template < class FT >
void
barycenterC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &p2z, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &p3z, const FT &w3,
             const FT &p4x, const FT &p4y, const FT &p4z, const FT &w4,
             FT &x, FT &y, FT &z)
{
   FT sum = w1 + w2 + w3 + w4;
   CGAL_kernel_assertion(sum != 0);
   x = (w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x) / sum;
   y = (w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y) / sum;
   z = (w1 * p1z + w2 * p2z + w3 * p3z + w4 * p4z) / sum;
}

template < class FT >
void
centroidC3( const FT &px, const FT &py, const FT &pz,
            const FT &qx, const FT &qy, const FT &qz,
            const FT &rx, const FT &ry, const FT &rz,
            const FT &sx, const FT &sy, const FT &sz,
            FT &x, FT &y, FT &z)
{
   x = (px + qx + rx + sx) / 4;
   y = (py + qy + ry + sy) / 4;
   z = (pz + qz + rz + sz) / 4;
}

template < class FT >
void
centroidC3( const FT &px, const FT &py, const FT &pz,
            const FT &qx, const FT &qy, const FT &qz,
            const FT &rx, const FT &ry, const FT &rz,
            FT &x, FT &y, FT &z)
{
   x = (px + qx + rx) / 3;
   y = (py + qy + ry) / 3;
   z = (pz + qz + rz) / 3;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
squared_radiusC3(const FT &px, const FT &py, const FT &pz,
                       const FT &qx, const FT &qy, const FT &qz,
                       const FT &rx, const FT &ry, const FT &rz,
                       const FT &sx, const FT &sy, const FT &sz)
{
  // Translate p to origin to simplify the expression.
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT qpz = qz-pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT rpz = rz-pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + CGAL_NTS square(rpz);
  FT spx = sx-px;
  FT spy = sy-py;
  FT spz = sz-pz;
  FT sp2 = CGAL_NTS square(spx) + CGAL_NTS square(spy) + CGAL_NTS square(spz);

  FT num_x = determinant(qpy,qpz,qp2,
                               rpy,rpz,rp2,
                               spy,spz,sp2);
  FT num_y = determinant(qpx,qpz,qp2,
                               rpx,rpz,rp2,
                               spx,spz,sp2);
  FT num_z = determinant(qpx,qpy,qp2,
                               rpx,rpy,rp2,
                               spx,spy,sp2);
  FT den   = determinant(qpx,qpy,qpz,
                               rpx,rpy,rpz,
                               spx,spy,spz);
  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );

  return (CGAL_NTS square(num_x) + CGAL_NTS square(num_y)
        + CGAL_NTS square(num_z)) / CGAL_NTS square(2 * den);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
squared_radiusC3(const FT &px, const FT &py, const FT &pz,
                       const FT &qx, const FT &qy, const FT &qz,
                       const FT &sx, const FT &sy, const FT &sz)
{
  // Translate s to origin to simplify the expression.
  FT psx = px-sx;
  FT psy = py-sy;
  FT psz = pz-sz;
  FT ps2 = CGAL_NTS square(psx) + CGAL_NTS square(psy) + CGAL_NTS square(psz);
  FT qsx = qx-sx;
  FT qsy = qy-sy;
  FT qsz = qz-sz;
  FT qs2 = CGAL_NTS square(qsx) + CGAL_NTS square(qsy) + CGAL_NTS square(qsz);
  FT rsx = psy*qsz-psz*qsy;
  FT rsy = psz*qsx-psx*qsz;
  FT rsz = psx*qsy-psy*qsx;

  FT num_x = ps2 * determinant(qsy,qsz,rsy,rsz)
	   - qs2 * determinant(psy,psz,rsy,rsz);
  FT num_y = ps2 * determinant(qsx,qsz,rsx,rsz)
	   - qs2 * determinant(psx,psz,rsx,rsz);
  FT num_z = ps2 * determinant(qsx,qsy,rsx,rsy)
	   - qs2 * determinant(psx,psy,rsx,rsy);

  FT den   = determinant(psx,psy,psz,
                               qsx,qsy,qsz,
                               rsx,rsy,rsz);

  CGAL_kernel_assertion( den != 0 );

  return (CGAL_NTS square(num_x) + CGAL_NTS square(num_y)
        + CGAL_NTS square(num_z)) / CGAL_NTS square(2 * den);
}

template <class FT> 
CGAL_KERNEL_MEDIUM_INLINE
void            
plane_from_pointsC3(const FT &px, const FT &py, const FT &pz,
                    const FT &qx, const FT &qy, const FT &qz,
                    const FT &rx, const FT &ry, const FT &rz,
		    FT &pa, FT &pb, FT &pc, FT &pd)
{
  FT rpx = px-rx;
  FT rpy = py-ry;
  FT rpz = pz-rz;
  FT rqx = qx-rx;
  FT rqy = qy-ry;
  FT rqz = qz-rz;
  // Cross product rp * rq
  pa = rpy*rqz - rqy*rpz;
  pb = rpz*rqx - rqz*rpx;
  pc = rpx*rqy - rqx*rpy;
  pd = - pa*rx - pb*ry - pc*rz;
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
void
plane_from_point_directionC3(const FT &px, const FT &py, const FT &pz,
                             const FT &dx, const FT &dy, const FT &dz,
                             FT &pa, FT &pb, FT &pc, FT &pd)
{
  // d is the normal direction
  pa = dx; pb = dy; pc = dz; pd = -dx*px - dy*py - dz*pz;
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
void
point_on_planeC3(const FT &pa, const FT &pb, const FT &pc, const FT &pd,
                 FT &x, FT &y, FT &z)
{
  x = y = z = 0;
  if (! CGAL_NTS is_zero(pa))      x = -pd/pa;
  else if (! CGAL_NTS is_zero(pb)) y = -pd/pb;
  else                             z = -pd/pc;
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
void
projection_planeC3(const FT &pa, const FT &pb, const FT &pc, const FT &pd,
                   const FT &px, const FT &py, const FT &pz,
                   FT &x, FT &y, FT &z)
{
  // the equation of the plane is Ax+By+Cz+D=0
  // the normal direction is (A,B,C)
  // the projected point is p-lambda(A,B,C) where
  // A(x-lambda A) + B(y-lambda B) + C(z-lambda C) + D = 0

  FT num = pa*px + pb*py + pc*pz + pd;
  FT den = pa*pa + pb*pb + pc*pc;
  FT lambda = num / den;

  x = px - lambda * pa;
  y = py - lambda * pb;
  z = pz - lambda * pc;
}

template < class FT >
CGAL_KERNEL_INLINE
FT
squared_distanceC3( const FT &px, const FT &py, const FT &pz,
                    const FT &qx, const FT &qy, const FT &qz)
{
  return CGAL_NTS square(px-qx) + CGAL_NTS square(py-qy) +
	 CGAL_NTS square(pz-qz);
}

template < class FT >
CGAL_KERNEL_INLINE
FT
squared_radiusC3( const FT &px, const FT &py, const FT &pz,
                  const FT &qx, const FT &qy, const FT &qz)
{
  return squared_distanceC3(px, py, pz, qx, qy, qz) / 4;
}

template < class FT >
CGAL_KERNEL_INLINE
FT
scaled_distance_to_directionC3(const FT &pa, const FT &pb, const FT &pc,
                               const FT &px, const FT &py, const FT &pz)
{
  return pa*px + pb*py + pc*pz;
}

template < class FT >
CGAL_KERNEL_INLINE
FT
scaled_distance_to_planeC3(
     const FT &pa, const FT &pb, const FT &pc, const FT &pd,
     const FT &px, const FT &py, const FT &pz)
{
  return pa*px + pb*py + pc*pz + pd;
}

template < class FT >
CGAL_KERNEL_INLINE
FT
scaled_distance_to_planeC3(
     const FT &ppx, const FT &ppy, const FT &ppz,
     const FT &pqx, const FT &pqy, const FT &pqz,
     const FT &prx, const FT &pry, const FT &prz,
     const FT &px,  const FT &py,  const FT &pz)
{
  return determinant(ppx-px,ppy-py,ppz-pz,
                           pqx-px,pqy-py,pqz-pz,
                           prx-px,pry-py,prz-pz);
}

template < class FT >
CGAL_KERNEL_INLINE
void
bisector_of_pointsC3(const FT &px, const FT &py, const FT &pz,
		     const FT &qx, const FT &qy, const FT &qz,
		     FT &a, FT &b, FT &c, FT &d)
{
  a = 2*(px - qx);
  b = 2*(py - qy);
  c = 2*(pz - qz);
  d = CGAL_NTS square(qx) + CGAL_NTS square(qy) + CGAL_NTS square(qz)
    - CGAL_NTS square(px) - CGAL_NTS square(py) - CGAL_NTS square(pz);
}

template < class FT >
void
bisector_of_planesC3(const FT &pa, const FT &pb, const FT &pc, const FT &pd,
		     const FT &qa, const FT &qb, const FT &qc, const FT &qd,
		     FT &a, FT &b, FT &c, FT &d)
{
  // We normalize the equations of the 2 planes, and we then add them.
  FT n1 = CGAL_NTS sqrt(CGAL_NTS square(pa) + CGAL_NTS square(pb) +
                        CGAL_NTS square(pc));
  FT n2 = CGAL_NTS sqrt(CGAL_NTS square(qa) + CGAL_NTS square(qb) +
                        CGAL_NTS square(qc));
  a = n2 * pa + n1 * qa;
  b = n2 * pb + n1 * qb;
  c = n2 * pc + n1 * qc;
  d = n2 * pd + n1 * qd;

  // Care must be taken for the case when this produces a degenerate line.
  if (a == 0 && b == 0 && c == 0) {
    a = n2 * pa - n1 * qa;
    b = n2 * pb - n1 * qb;
    c = n2 * pc - n1 * qc;
    d = n2 * pd - n1 * qd;
  }
}

template < class FT >
FT
squared_areaC3(const FT &px, const FT &py, const FT &pz,
	       const FT &qx, const FT &qy, const FT &qz,
	       const FT &rx, const FT &ry, const FT &rz)
{
    // Compute vectors pq and pr, then the cross product,
    // then 1/4 of its squared length.
    FT dqx = qx-px;
    FT dqy = qy-py;
    FT dqz = qz-pz;
    FT drx = rx-px;
    FT dry = ry-py;
    FT drz = rz-pz;

    FT vx = dqy*drz-dqz*dry;
    FT vy = dqz*drx-dqx*drz;
    FT vz = dqx*dry-dqy*drx;

    return (CGAL_NTS square(vx) + CGAL_NTS square(vy) + CGAL_NTS square(vz))/4;
}

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
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz) - qw + pw;
  FT rpx = rx - px;
  FT rpy = ry - py;
  FT rpz = rz - pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) +
           CGAL_NTS square(rpz) - rw + pw;
  FT spx = sx - px;
  FT spy = sy - py;
  FT spz = sz - pz;
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

template <class FT>
void
determinants_for_circumcenterC3(const FT &px, const FT &py, const FT &pz,
                                const FT &qx, const FT &qy, const FT &qz,
                                const FT &rx, const FT &ry, const FT &rz,
                                const FT &sx, const FT &sy, const FT &sz,
                                FT &num_x,  FT &num_y, FT &num_z, FT& den)
{
  // translate origin to p
  // and compute determinants for weighted_circumcenter and
  // circumradius
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT rpx = rx - px;
  FT rpy = ry - py;
  FT rpz = rz - pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) +
           CGAL_NTS square(rpz);
  FT spx = sx - px;
  FT spy = sy - py;
  FT spz = sz - pz;
  FT sp2 = CGAL_NTS square(spx) + CGAL_NTS square(spy) +
           CGAL_NTS square(spz);

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
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                        const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                        const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                        const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                        FT &x, FT &y, FT &z)
{
  // this function computes the weighted circumcenter point only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
                                           qx, qy, qz, qw,
                                           rx, ry, rz, rw,
                                           sx, sy, sz, sw,
                                           num_x,  num_y, num_z,den);

  CGAL_assertion( ! CGAL_NTS is_zero(den) );
  FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
}

template < class FT>
void
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                        const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                        const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                        const FT &sx, const FT &sy, const FT &sz, const FT &sw,
                        FT &x, FT &y, FT &z, FT &w)
{
  // this function computes the weighted circumcenter point
  // and the squared weighted circumradius

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
                                           qx, qy, qz, qw,
                                           rx, ry, rz, rw,
                                           sx, sy, sz, sw,
                                           num_x,  num_y, num_z, den);

  CGAL_assertion( ! CGAL_NTS is_zero(den) );
  FT inv = FT(1)/(FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;

  w = (CGAL_NTS square(num_x) + CGAL_NTS square(num_y) + CGAL_NTS square(num_z))
      * CGAL_NTS square(inv) - pw;
}

template< class FT >
FT
squared_radius_orthogonal_sphereC3(
    const FT &px, const FT &py, const FT &pz, const FT &pw,
    const FT &qx, const FT &qy, const FT &qz, const FT &qw,
    const FT &rx, const FT &ry, const FT &rz, const FT &rw,
    const FT &sx, const FT &sy, const FT &sz, const FT &sw)
{
  // this function computes the squared weighted circumradius only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
                                           qx, qy, qz, qw,
                                           rx, ry, rz, rw,
                                           sx, sy, sz, sw,
                                           num_x, num_y, num_z,den);

  CGAL_assertion( ! CGAL_NTS is_zero(den) );
  FT inv = FT(1)/(FT(2) * den);

  return (CGAL_NTS square(num_x) + CGAL_NTS square(num_y) + CGAL_NTS square(num_z))
         * CGAL_NTS square(inv) - pw;
}

template <class FT>
void
determinants_for_weighted_circumcenterC3(
    const FT &px, const FT &py, const FT &pz, const FT &pw,
    const FT &qx, const FT &qy, const FT &qz, const FT &qw,
    const FT &rx, const FT &ry, const FT &rz, const FT &rw,
    FT &num_x,  FT &num_y, FT &num_z, FT &den)
{
  // translate origin to p and compute determinants for weighted_circumcenter
  // and circumradius

  // Translate s to origin to simplify the expression.
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz) - qw + pw;
  FT rpx = rx - px;
  FT rpy = ry - py;
  FT rpz = rz - pz;
  FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) +
           CGAL_NTS square(rpz) - rw + pw;

  FT sx = qpy*rpz - qpz*rpy;
  FT sy = qpz*rpx - qpx*rpz;
  FT sz = qpx*rpy - qpy*rpx;

  // The following determinants can be developped and simplified.
//
//  FT num_x = determinant(qpy,qpz,qp2,
//                         rpy,rpz,rp2,
//                         sy,sz,FT(0));
//  FT num_y = determinant(qpx,qpz,qp2,
//                         rpx,rpz,rp2,
//                         sx,sz,FT(0));
//  FT num_z = determinant(qpx,qpy,qp2,
//                         rpx,rpy,rp2,
//                         sx,sy,FT(0));

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
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                        const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                        const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                        FT &x, FT &y, FT &z)
{
  // this function computes the weighted circumcenter point only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
                                           qx, qy, qz, qw,
                                           rx, ry, rz, rw,
                                           num_x,  num_y, num_z, den);

  CGAL_assertion( den != FT(0) );
  FT inv = FT(1) / (FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;
}

template < class FT >
void
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                        const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                        const FT &rx, const FT &ry, const FT &rz, const FT &rw,
                        FT &x, FT &y, FT &z, FT &w)
{
  // this function computes the weighted circumcenter and
  // the weighted squared circumradius

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
                                           qx, qy, qz, qw,
                                           rx, ry, rz, rw,
                                           num_x,  num_y, num_z, den);

  CGAL_assertion( den != FT(0) );
  FT inv = FT(1) / (FT(2) * den);

  x = px + num_x*inv;
  y = py - num_y*inv;
  z = pz + num_z*inv;

  w = (CGAL_NTS square(num_x) + CGAL_NTS square(num_y) + CGAL_NTS square(num_z))
      *CGAL_NTS square(inv) - pw;
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
    const FT &px, const FT &py, const FT &pz, const FT &pw,
    const FT &qx, const FT &qy, const FT &qz, const FT &qw,
    const FT &rx, const FT &ry, const FT &rz, const FT &rw)
{
  // this function computes the weighted squared circumradius only

  // Translate p to origin and compute determinants
  FT num_x, num_y, num_z, den;
  determinants_for_weighted_circumcenterC3(px, py, pz, pw,
                                           qx, qy, qz, qw,
                                           rx, ry, rz, rw,
                                           num_x,  num_y, num_z, den);

  CGAL_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

  return (CGAL_NTS square(num_x) + CGAL_NTS square(num_y) + CGAL_NTS square(num_z))
         * CGAL_NTS square(inv) - pw;
}

template < class FT >
void
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                        const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                        FT &x, FT &y, FT &z)
{
// this function computes the weighted circumcenter point only
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT inv = FT(1) / (FT(2) * qp2);
  FT alpha = 1 / FT(2) + (pw-qw) * inv;

  x = px + alpha * qpx;
  y = py + alpha * qpy;
  z = pz + alpha * qpz;
}

template < class FT >
void
weighted_circumcenterC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                        const FT &qx, const FT &qy, const FT &qz, const FT &qw,
                        FT &x, FT &y, FT &z, FT &w)
{
 // this function computes the weighted circumcenter point and
  // the weighted circumradius
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT inv = FT(1) / (FT(2) * qp2);
  FT alpha = 1 / FT(2) + (pw-qw) * inv;

  x = px + alpha * qpx;
  y = py + alpha * qpy;
  z = pz + alpha * qpz;

  w = CGAL_NTS square(alpha) * qp2 - pw;
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_sphereC3(
  const FT &px, const FT &py, const FT &pz, const FT  &pw,
  const FT &qx, const FT &qy, const FT &qz, const FT  &qw)
{
  // this function computes the weighted circumradius only
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) +
           CGAL_NTS square(qpz);
  FT inv = FT(1) / (FT(2) * qp2);
  FT alpha = 1 / FT(2) + (pw-qw) * inv;

  return  CGAL_NTS square(alpha)*qp2 - pw;
}

template< class FT >
FT
power_productC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
                const FT &qx, const FT &qy, const FT &qz, const FT &qw)
{
  // computes the power product of two weighted points
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qpz = qz - pz;
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

  a = RT(1)*determinant(dqy, dqz, dry, drz);
  b = - RT(1)*determinant(dqx, dqz, drx, drz);
  c = RT(1)*determinant(dqx, dqy, drx, dry);
}

// function used in critical_squared_radiusC3
// power ( t, tw) with respect to
// circle orthogonal (p,pw), (q,qw), (r,rw), (s,sw)
template < class FT>
FT
power_to_orthogonal_sphereC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
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
power_distance_to_power_sphereC3(const FT &px, const FT &py, const FT &pz, const FT &pw,
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

  return -ff0/(ff1 - ff0);
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
//           + FT2*Lambda*Mu*dpdq - rw);

} //namespace CGAL

#endif // CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
