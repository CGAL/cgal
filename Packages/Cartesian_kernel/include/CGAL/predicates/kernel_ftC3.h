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
// file          : include/CGAL/predicates/kernel_ftC3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_PREDICATES_KERNEL_FTC3_H
#define CGAL_PREDICATES_KERNEL_FTC3_H

#include <CGAL/predicates/sign_of_determinant.h>
#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_lexicographically_xyzC3(const FT &px, const FT &py, const FT &pz,
                                const FT &qx, const FT &qy, const FT &qz)
{
  Comparison_result c = CGAL_NTS compare(px, qx);
  if (c != EQUAL) return c;
  c = CGAL_NTS compare(py, qy);
  if (c != EQUAL) return c;
  return CGAL_NTS compare(pz, qz);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
strict_dominanceC3(const FT &px, const FT &py, const FT &pz,
		   const FT &qx, const FT &qy, const FT &qz)
{
  return CGAL_NTS compare(px, qx) == LARGER &&
	 CGAL_NTS compare(py, qy) == LARGER &&
	 CGAL_NTS compare(pz, qz) == LARGER;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
dominanceC3(const FT &px, const FT &py, const FT &pz,
	    const FT &qx, const FT &qy, const FT &qz)
{
  return CGAL_NTS compare(px, qx) != SMALLER && 
	 CGAL_NTS compare(py, qy) != SMALLER &&
	 CGAL_NTS compare(pz, qz) != SMALLER;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
collinearC3(const FT &px, const FT &py, const FT &pz,
            const FT &qx, const FT &qy, const FT &qz,
            const FT &rx, const FT &ry, const FT &rz)
{
  FT dpx = px-rx;
  FT dqx = qx-rx;
  FT dpy = py-ry;
  FT dqy = qy-ry;
  if (sign_of_determinant2x2(dpx, dqx, dpy, dqy) != ZERO)
      return false;
  FT dpz = pz-rz;
  FT dqz = qz-rz;
  return sign_of_determinant2x2(dpx, dqx, dpz, dqz) == ZERO
      && sign_of_determinant2x2(dpy, dqy, dpz, dqz) == ZERO;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Orientation
orientationC3(const FT &px, const FT &py, const FT &pz,
              const FT &qx, const FT &qy, const FT &qz,
              const FT &rx, const FT &ry, const FT &rz,
              const FT &sx, const FT &sy, const FT &sz)
{
  return Orientation(sign_of_determinant3x3(qx-px,rx-px,sx-px,
                                            qy-py,ry-py,sy-py,
                                            qz-pz,rz-pz,sz-pz));
}

template < class FT >
inline
Angle
angleC3(const FT &px, const FT &py, const FT &pz,
        const FT &qx, const FT &qy, const FT &qz,
        const FT &rx, const FT &ry, const FT &rz)
{
  return (Angle) CGAL_NTS sign ((px-qx)*(rx-qx)+
	                        (py-qy)*(ry-qy)+
				(pz-qz)*(rz-qz));
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
Orientation
coplanar_orientationC3(const FT &px, const FT &py, const FT &pz,
                       const FT &qx, const FT &qy, const FT &qz,
                       const FT &rx, const FT &ry, const FT &rz,
                       const FT &sx, const FT &sy, const FT &sz)
{
  Orientation oxy_pqr = orientationC2(px,py,qx,qy,rx,ry);
  if (oxy_pqr != COLLINEAR)
      return Orientation( oxy_pqr * orientationC2(px,py,qx,qy,sx,sy));

  Orientation oyz_pqr = orientationC2(py,pz,qy,qz,ry,rz);
  if (oyz_pqr != COLLINEAR)
      return Orientation( oyz_pqr * orientationC2(py,pz,qy,qz,sy,sz));

  Orientation oxz_pqr = orientationC2(px,pz,qx,qz,rx,rz);
  CGAL_kernel_assertion(oxz_pqr != COLLINEAR);
  return Orientation( oxz_pqr * orientationC2(px,pz,qx,qz,sx,sz));
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
Orientation
coplanar_orientationC3(const FT &px, const FT &py, const FT &pz,
                       const FT &qx, const FT &qy, const FT &qz,
                       const FT &rx, const FT &ry, const FT &rz)
{
  Orientation oxy_pqr = orientationC2(px,py,qx,qy,rx,ry);
  if (oxy_pqr != COLLINEAR)
      return oxy_pqr;

  Orientation oyz_pqr = orientationC2(py,pz,qy,qz,ry,rz);
  if (oyz_pqr != COLLINEAR)
      return oyz_pqr;

  return orientationC2(px,pz,qx,qz,rx,rz);
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
coplanar_side_of_bounded_circleC3(const FT &px, const FT &py, const FT &pz,
                                  const FT &qx, const FT &qy, const FT &qz,
                                  const FT &rx, const FT &ry, const FT &rz,
                                  const FT &tx, const FT &ty, const FT &tz)
{
  // The approach is to compute side_of_bounded_sphere(p,q,r,t+v,t),
  // with v = pq ^ pr.
  // Note : since the circle defines the orientation of the plane, it can not
  // be considered oriented.
  FT ptx = px - tx;
  FT pty = py - ty;
  FT ptz = pz - tz;
  FT pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty) + CGAL_NTS square(ptz);
  FT qtx = qx - tx;
  FT qty = qy - ty;
  FT qtz = qz - tz;
  FT qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty) + CGAL_NTS square(qtz);
  FT rtx = rx - tx;
  FT rty = ry - ty;
  FT rtz = rz - tz;
  FT rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty) + CGAL_NTS square(rtz);
  FT pqx = qx - px;
  FT pqy = qy - py;
  FT pqz = qz - pz;
  FT prx = rx - px;
  FT pry = ry - py;
  FT prz = rz - pz;
  FT vx = pqy*prz - pqz*pry;
  FT vy = pqz*prx - pqx*prz;
  FT vz = pqx*pry - pqy*prx;
  FT v2 = CGAL_NTS square(vx) + CGAL_NTS square(vy) + CGAL_NTS square(vz);
  return Bounded_side(sign_of_determinant4x4(ptx,pty,ptz,pt2,
                                             rtx,rty,rtz,rt2,
                                             qtx,qty,qtz,qt2,
                                             vx,vy,vz,v2));
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear_are_ordered_along_lineC3(
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz,
     const FT &rx, const FT &ry, const FT &rz)
{
  if (px < qx) return !(rx < qx);
  if (qx < px) return !(qx < rx);
  if (py < qy) return !(ry < qy);
  if (qy < py) return !(qy < ry);
  if (pz < qz) return !(rz < qz);
  if (qz < pz) return !(qz < rz);
  return true; // p==q
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear_are_strictly_ordered_along_lineC3(
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz,
     const FT &rx, const FT &ry, const FT &rz)
{
  if (px < qx) return (qx < rx);
  if (qx < px) return (rx < qx);
  if (py < qy) return (qy < ry);
  if (qy < py) return (ry < qy);
  if (pz < qz) return (qz < rz);
  if (qz < pz) return (rz < qz);
  return false; // p==q
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
equal_directionC3(const FT &dx1, const FT &dy1, const FT &dz1,
                  const FT &dx2, const FT &dy2, const FT &dz2)
{
  return sign_of_determinant2x2(dx1, dy1, dx2, dy2) == ZERO
      && sign_of_determinant2x2(dx1, dz1, dx2, dz2) == ZERO
      && sign_of_determinant2x2(dy1, dz1, dy2, dz2) == ZERO
      && CGAL_NTS sign(dx1) == CGAL_NTS sign(dx2)
      && CGAL_NTS sign(dy1) == CGAL_NTS sign(dy2)
      && CGAL_NTS sign(dz1) == CGAL_NTS sign(dz2);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
equal_planeC3(const FT &ha, const FT &hb, const FT &hc, const FT &hd,
              const FT &pa, const FT &pb, const FT &pc, const FT &pd)
{
    if (!equal_directionC3(ha, hb, hc, pa, pb, pc))
	return false; // Not parallel.

    CGAL::Sign s1a = CGAL_NTS sign(ha);
    if (s1a != ZERO)
        return s1a == CGAL_NTS sign(pa)
            && sign_of_determinant2x2(pa, pd, ha, hd) == ZERO;
    CGAL::Sign s1b = CGAL_NTS sign(hb);
    if (s1b != ZERO)
        return s1b == CGAL_NTS sign(pb)
            && sign_of_determinant2x2(pb, pd, hb, hd) == ZERO;
    return CGAL_NTS sign(pc) == CGAL_NTS sign(hc)
        && sign_of_determinant2x2(pc, pd, hc, hd) == ZERO;
}

template <class FT >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
side_of_oriented_planeC3(const FT &a,  const FT &b,  const FT &c, const FT &d,
                         const FT &px, const FT &py, const FT &pz)
{
  return Oriented_side(CGAL_NTS sign(a*px + b*py + c*pz +d));
}

template <class FT >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
side_of_oriented_sphereC3(const FT &px, const FT &py, const FT &pz,
                          const FT &qx, const FT &qy, const FT &qz,
                          const FT &rx, const FT &ry, const FT &rz,
                          const FT &sx, const FT &sy, const FT &sz,
                          const FT &tx, const FT &ty, const FT &tz)
{
  FT ptx = px - tx;
  FT pty = py - ty;
  FT ptz = pz - tz;
  FT pt2 = CGAL_NTS square(ptx) + CGAL_NTS square(pty) + CGAL_NTS square(ptz);
  FT qtx = qx - tx;
  FT qty = qy - ty;
  FT qtz = qz - tz;
  FT qt2 = CGAL_NTS square(qtx) + CGAL_NTS square(qty) + CGAL_NTS square(qtz);
  FT rtx = rx - tx;
  FT rty = ry - ty;
  FT rtz = rz - tz;
  FT rt2 = CGAL_NTS square(rtx) + CGAL_NTS square(rty) + CGAL_NTS square(rtz);
  FT stx = sx - tx;
  FT sty = sy - ty;
  FT stz = sz - tz;
  FT st2 = CGAL_NTS square(stx) + CGAL_NTS square(sty) + CGAL_NTS square(stz);
  return Oriented_side(sign_of_determinant4x4(ptx,pty,ptz,pt2,
                                              rtx,rty,rtz,rt2,
                                              qtx,qty,qtz,qt2,
                                              stx,sty,stz,st2));
}

template <class FT >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
side_of_bounded_sphereC3(const FT &px, const FT &py, const FT &pz,
                         const FT &qx, const FT &qy, const FT &qz,
                         const FT &rx, const FT &ry, const FT &rz,
                         const FT &sx, const FT &sy, const FT &sz,
                         const FT &tx, const FT &ty, const FT &tz)
{
  Oriented_side s = side_of_oriented_sphereC3(px, py, pz,
                                              qx, qy, qz,
                                              rx, ry, rz,
                                              sx, sy, sz,
                                              tx, ty, tz);
  Orientation o = orientationC3(px, py, pz,
                                qx, qy, qz,
                                rx, ry, rz,
                                sx, sy, sz);
  return Bounded_side(s * o);
}

template <class FT >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
side_of_bounded_sphereC3(const FT &px, const FT &py, const FT &pz,
                         const FT &qx, const FT &qy, const FT &qz,
                         const FT &tx, const FT &ty, const FT &tz)
{
  // Returns whether T lies inside or outside the sphere which diameter is PQ.
  return Bounded_side( CGAL_NTS sign((tx-px)*(qx-tx)
	                           + (ty-py)*(qy-ty)
	                           + (tz-pz)*(qz-tz)) );
}

template < class FT >
CGAL_KERNEL_INLINE
Comparison_result
cmp_dist_to_pointC3(const FT &px, const FT &py, const FT &pz,
                    const FT &qx, const FT &qy, const FT &qz,
                    const FT &rx, const FT &ry, const FT &rz)
{
  return CGAL_NTS compare(squared_distanceC3(px,py,pz,qx,qy,qz),
                          squared_distanceC3(px,py,pz,rx,ry,rz));
}

// Because of the way the filtered predicates generator script works,
// cmp_dist_to_pointC3() must be defined _before_ ths following one.
template <class FT >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
side_of_bounded_sphereC3(const FT &px, const FT &py, const FT &pz,
                         const FT &qx, const FT &qy, const FT &qz,
                         const FT &sx, const FT &sy, const FT &sz,
                         const FT &tx, const FT &ty, const FT &tz)
{
  // Returns whether T lies inside or outside the sphere which equatorial
  // circle is PQR.

  // This code is inspired by the one of circumcenterC3(3 points).

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
  FT tsx = tx-sx;
  FT tsy = ty-sy;
  FT tsz = tz-sz;

  FT num_x = ps2 * det2x2_by_formula(qsy,qsz,rsy,rsz)
	   - qs2 * det2x2_by_formula(psy,psz,rsy,rsz);
  FT num_y = ps2 * det2x2_by_formula(qsx,qsz,rsx,rsz)
	   - qs2 * det2x2_by_formula(psx,psz,rsx,rsz);
  FT num_z = ps2 * det2x2_by_formula(qsx,qsy,rsx,rsy)
	   - qs2 * det2x2_by_formula(psx,psy,rsx,rsy);

  FT den   = det3x3_by_formula(psx,psy,psz,
                               qsx,qsy,qsz,
                               rsx,rsy,rsz);

  FT den2 = FT(2) * den;

  // The following could be simplified a bit.
  return Bounded_side(cmp_dist_to_pointC3(num_x,    - num_y,  num_z,
	                                  psx*den2, psy*den2, psz*den2,
	                                  tsx*den2, tsy*den2, tsz*den2) );
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_dist_to_pointC3(const FT &px, const FT &py, const FT &pz,
                           const FT &qx, const FT &qy, const FT &qz,
                           const FT &rx, const FT &ry, const FT &rz)
{
  return cmp_dist_to_pointC3(px,py,pz,qx,qy,qz,rx,ry,rz) == LARGER;
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_dist_to_pointC3(const FT &px, const FT &py, const FT &pz,
                            const FT &qx, const FT &qy, const FT &qz,
                            const FT &rx, const FT &ry, const FT &rz)
{
  return cmp_dist_to_pointC3(px,py,pz,qx,qy,qz,rx,ry,rz) == SMALLER;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
cmp_signed_dist_to_directionC3(
     const FT &pa, const FT &pb, const FT &pc,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return CGAL_NTS compare(scaled_distance_to_directionC3(pa,pb,pc,px,py,pz),
                          scaled_distance_to_directionC3(pa,pb,pc,qx,qy,qz));
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_signed_dist_to_directionC3(
     const FT &pa, const FT &pb, const FT &pc,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return cmp_signed_dist_to_directionC3(pa,pb,pc,px,py,pz,qx,qy,qz) == LARGER;
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_signed_dist_to_directionC3(
     const FT &pa, const FT &pb, const FT &pc,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return cmp_signed_dist_to_directionC3(pa,pb,pc,px,py,pz,qx,qy,qz) == SMALLER;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
cmp_signed_dist_to_planeC3(
     const FT &ppx, const FT &ppy, const FT &ppz,
     const FT &pqx, const FT &pqy, const FT &pqz,
     const FT &prx, const FT &pry, const FT &prz,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return Comparison_result(sign_of_determinant3x3(
	      pqx-ppx, pqy-ppy, pqz-ppz,
	      prx-ppx, pry-ppy, prz-ppz,
	      qx-px,   qy-py,   qz-pz));
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_signed_dist_to_planeC3(
     const FT &ppx, const FT &ppy, const FT &ppz,
     const FT &pqx, const FT &pqy, const FT &pqz,
     const FT &prx, const FT &pry, const FT &prz,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
    return cmp_signed_dist_to_planeC3(ppx, ppy, ppz, pqx, pqy, pqz,
	    prx, pry, prz, px, py, pz, qx, qy, qz) == LARGER;
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_signed_dist_to_planeC3(
     const FT &ppx, const FT &ppy, const FT &ppz,
     const FT &pqx, const FT &pqy, const FT &pqz,
     const FT &prx, const FT &pry, const FT &prz,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
    return cmp_signed_dist_to_planeC3(ppx, ppy, ppz, pqx, pqy, pqz,
	    prx, pry, prz, px, py, pz, qx, qy, qz) == SMALLER;
}

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/predicates/kernel_ftC3.h>
#endif

#endif // CGAL_PREDICATES_KERNEL_FTC3_H
