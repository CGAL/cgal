// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/constructions/kernel_ftC3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
#define CGAL_CONSTRUCTIONS_KERNEL_FTC3_H

#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_INLINE
void
midpointC3( const FT &px, const FT &py, const FT &pz,
            const FT &qx, const FT &qy, const FT &qz,
            FT &x, FT &y, FT &z)
{
  FT half = FT(1) / FT(2);
  x = (px+qx) * half;
  y = (py+qy) * half;
  z = (pz+qz) * half;
}

template < class FT >
void
circumcenterC3( const FT &px, const FT &py, const FT &pz,
                const FT &qx, const FT &qy, const FT &qz,
                const FT &rx, const FT &ry, const FT &rz,
                const FT &sx, const FT &sy, const FT &sz,
                FT &x, FT &y, FT &z)
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
  FT rsx = rx-sx;
  FT rsy = ry-sy;
  FT rsz = rz-sz;
  FT rs2 = CGAL_NTS square(rsx) + CGAL_NTS square(rsy) + CGAL_NTS square(rsz);

  FT num_x = det3x3_by_formula(psy,psz,ps2,
                               qsy,qsz,qs2,
                               rsy,rsz,rs2);
  FT num_y = det3x3_by_formula(psx,psz,ps2,
                               qsx,qsz,qs2,
                               rsx,rsz,rs2);
  FT num_z = det3x3_by_formula(psx,psy,ps2,
                               qsx,qsy,qs2,
                               rsx,rsy,rs2);
  FT den   = det3x3_by_formula(psx,psy,psz,
                               qsx,qsy,qsz,
                               rsx,rsy,rsz);
  CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );
  FT inv = FT(1)/(FT(2) * den);

  x = sx + num_x*inv;
  y = sy - num_y*inv;
  z = sz + num_z*inv;
}

template < class FT >
void
circumcenterC3( const FT &px, const FT &py, const FT &pz,
                const FT &qx, const FT &qy, const FT &qz,
                const FT &sx, const FT &sy, const FT &sz,
                FT &x, FT &y, FT &z)
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

  // The following determinants can be developped and simplified.
  //
  // FT num_x = det3x3_by_formula(psy,psz,ps2,
  //                              qsy,qsz,qs2,
  //                              rsy,rsz,FT(0));
  // FT num_y = det3x3_by_formula(psx,psz,ps2,
  //                              qsx,qsz,qs2,
  //                              rsx,rsz,FT(0));
  // FT num_z = det3x3_by_formula(psx,psy,ps2,
  //                              qsx,qsy,qs2,
  //                              rsx,rsy,FT(0));

  FT num_x = ps2 * det2x2_by_formula(qsy,qsz,rsy,rsz)
	   - qs2 * det2x2_by_formula(psy,psz,rsy,rsz);
  FT num_y = ps2 * det2x2_by_formula(qsx,qsz,rsx,rsz)
	   - qs2 * det2x2_by_formula(psx,psz,rsx,rsz);
  FT num_z = ps2 * det2x2_by_formula(qsx,qsy,rsx,rsy)
	   - qs2 * det2x2_by_formula(psx,psy,rsx,rsy);

  FT den   = det3x3_by_formula(psx,psy,psz,
                               qsx,qsy,qsz,
                               rsx,rsy,rsz);

  CGAL_kernel_assertion( den != FT(0) );
  FT inv = FT(1)/(FT(2) * den);

  x = sx + num_x*inv;
  y = sy - num_y*inv;
  z = sz + num_z*inv;
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
void
point_on_lineC3(const FT &lpx, const FT &lpy, const FT &lpz,
                const FT &ldx, const FT &ldy, const FT &ldz,
		int i,
                FT &x, FT &y, FT &z)
{
  x = lpx + ldx*FT(i);
  y = lpy + ldy*FT(i);
  z = lpz + ldz*FT(i);
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
void
projection_lineC3(const FT &px, const FT &py, const FT &pz,
                  const FT &lpx, const FT &lpy, const FT &lpz,
                  const FT &ldx, const FT &ldy, const FT &ldz,
                  FT &x, FT &y, FT &z)
{
  // projects p on the line l
  FT dpx = px-lpx;
  FT dpy = py-lpy;
  FT dpz = pz-lpz;
  FT lambda = (ldx*dpx+ldy*dpy+ldz*dpz) / (ldx*ldx+ldy*ldy+ldz*ldz);
  x = lpx + lambda * ldx;
  y = lpy + lambda * ldy;
  z = lpz + lambda * ldz;
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
  x = y = z = FT(0);
  if (! CGAL_NTS is_zero(pa))      x = -pd/pa;
  else if (! CGAL_NTS is_zero(pb)) y = -pd/pb;
  else                  z = -pd/pc;
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
  return det3x3_by_formula(ppx-px,ppy-py,ppz-pz,
                           pqx-px,pqy-py,pqz-pz,
                           prx-px,pry-py,prz-pz);
}

CGAL_END_NAMESPACE

#endif // CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
