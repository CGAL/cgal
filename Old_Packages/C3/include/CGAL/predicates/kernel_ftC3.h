// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//

// release       :
// release_date  :
//
// file          : include/CGAL/predicates/kernel_ftC3.h
// source        : include/CGAL/predicates/kernel_ftC3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_PREDICATES_KERNEL_FTC3_H
#define CGAL_PREDICATES_KERNEL_FTC3_H

#ifndef CGAL_DETERMINANT_H
#include <CGAL/determinant.h>
#endif // CGAL_DETERMINANT_H
#ifndef CGAL_PREDICATES_KERNEL_FTC2_H
#include <CGAL/predicates/kernel_ftC2.h>
#endif // CGAL_PREDICATES_KERNEL_FTC2_H
#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
#include <CGAL/constructions/kernel_ftC3.h>
#endif // CGAL_CONSTRUCTIONS_KERNEL_FTC3_H

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_lexicographically_xyzC3(const FT &x1, const FT &y1, const FT &z1,
                                const FT &x2, const FT &y2, const FT &z2)
{
  if (x1 < x2) return SMALLER;
  if (x2 < x1) return LARGER;
  if (y1 < y2) return SMALLER;
  if (y2 < y1) return LARGER;
  if (z1 < z2) return SMALLER;
  if (z2 < z1) return LARGER;
  return EQUAL;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
collinearC3(const FT &px, const FT &py, const FT &pz,
            const FT &qx, const FT &qy, const FT &qz,
            const FT &rx, const FT &ry, const FT &rz)
{
  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dpz = pz-rz;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dqz = qz-rz;
  return (sign_of_determinant2x2(dpx,dqx,dpy,dqy) == ZERO)
      && (sign_of_determinant2x2(dpx,dqx,dpz,dqz) == ZERO)
      && (sign_of_determinant2x2(dpy,dqy,dpz,dqz) == ZERO);
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
CGAL_KERNEL_MEDIUM_INLINE
Orientation
coplanar_orientationC3(const FT &qx, const FT &qy, const FT &qz,
                       const FT &rx, const FT &ry, const FT &rz,
                       const FT &sx, const FT &sy, const FT &sz,
                       const FT &px, const FT &py, const FT &pz)
{
  Orientation oxy_qrs = orientationC2(qx,qy,rx,ry,sx,sy);
  if (oxy_qrs != COLLINEAR)
      return Orientation( oxy_qrs * orientationC2(qx,qy,rx,ry,px,py));
  Orientation oyz_qrs = orientationC2(qy,qz,ry,rz,sy,sz);
  if (oyz_qrs != COLLINEAR)
      return Orientation( oyz_qrs * orientationC2(qy,qz,ry,rz,py,pz));
  Orientation oxz_qrs = orientationC2(qx,qz,rx,rz,sx,sz);
  assert(oxz_qrs != COLLINEAR);
  return Orientation( oxz_qrs * orientationC2(qx,qz,rx,rz,px,pz));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear_are_ordered_along_lineC3(
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz,
     const FT &rx, const FT &ry, const FT &rz)
{
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
  if (px < qx) return !(rx < qx);
  if (qx < px) return !(qx < rx);
  if (py < qy) return !(ry < qy);
  if (qy < py) return !(qy < ry);
  if (pz < qz) return !(rz < qz);
  if (qz < pz) return !(qz < rz);
  return true; // p==q
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
collinear_are_strictly_ordered_along_lineC3(
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz,
     const FT &rx, const FT &ry, const FT &rz)
{
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
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
  return ( dx1*dy2 == dy1*dx2 )
       &&( dx1*dz2 == dz1*dx2 )
       &&( dy1*dz2 == dz1*dy2 )
       &&( CGAL::sign(dx1) == CGAL::sign(dx2) )
       &&( CGAL::sign(dy1) == CGAL::sign(dy2) )
       &&( CGAL::sign(dz1) == CGAL::sign(dz2) );
}

template <class FT >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
side_of_oriented_planeC3(const FT &a,const FT &b,const FT &c,const FT &d,
                         const FT &px, const FT &py, const FT &pz)
{
  return Oriented_side(CGAL::sign(a*px + b*py + c*pz +d));
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
  FT pt2 = square(ptx) + square(pty) + square(ptz);
  FT qtx = qx - tx;
  FT qty = qy - ty;
  FT qtz = qz - tz;
  FT qt2 = square(qtx) + square(qty) + square(qtz);
  FT rtx = rx - tx;
  FT rty = ry - ty;
  FT rtz = rz - tz;
  FT rt2 = square(rtx) + square(rty) + square(rtz);
  FT stx = sx - tx;
  FT sty = sy - ty;
  FT stz = sz - tz;
  FT st2 = square(stx) + square(sty) + square(stz);
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


template < class FT >
CGAL_KERNEL_INLINE
Comparison_result
cmp_dist_to_pointC3(const FT &px, const FT &py, const FT &pz,
                    const FT &qx, const FT &qy, const FT &qz,
                    const FT &rx, const FT &ry, const FT &rz)
{
  return CGAL::compare(squared_distanceC3(px,py,pz,qx,qy,qz),
                       squared_distanceC3(px,py,pz,rx,ry,rz));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_dist_to_pointC3(const FT &px, const FT &py, const FT &pz,
                           const FT &qx, const FT &qy, const FT &qz,
                           const FT &rx, const FT &ry, const FT &rz)
{
  return squared_distanceC3(px,py,pz,rx,ry,rz)
       < squared_distanceC3(px,py,pz,qx,qy,qz);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_dist_to_pointC3(const FT &px, const FT &py, const FT &pz,
                            const FT &qx, const FT &qy, const FT &qz,
                            const FT &rx, const FT &ry, const FT &rz)
{
  return squared_distanceC3(px,py,pz,qx,qy,qz)
       < squared_distanceC3(px,py,pz,rx,ry,rz);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
cmp_signed_dist_to_planeC3(
     const FT &pa, const FT &pb, const FT &pc, const FT &pd,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return CGAL::compare(scaled_distance_to_planeC3(pa,pb,pc,pd,px,py,pz),
                       scaled_distance_to_planeC3(pa,pb,pc,pd,qx,qy,qz));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_signed_dist_to_planeC3(
     const FT &pa, const FT &pb, const FT &pc, const FT &pd,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return scaled_distance_to_planeC3(pa,pb,pc,pd,qx,qy,qz)
       < scaled_distance_to_planeC3(pa,pb,pc,pd,px,py,pz);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_signed_dist_to_planeC3(
     const FT &pa, const FT &pb, const FT &pc, const FT &pd,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return scaled_distance_to_planeC3(pa,pb,pc,pd,px,py,pz)
       < scaled_distance_to_planeC3(pa,pb,pc,pd,qx,qy,qz);
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
  return CGAL::compare(
           scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                      prx,pry,prz,psx,psy,psz,
                                      px,py,pz),
           scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                      prx,pry,prz,psx,psy,psz,
                                      qx,qy,qz) );
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
has_larger_signed_dist_to_planeC3(
     const FT &ppx, const FT &ppy, const FT &ppz,
     const FT &pqx, const FT &pqy, const FT &pqz,
     const FT &prx, const FT &pry, const FT &prz,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                    prx,pry,prz,psx,psy,psz,
                                    qx,qy,qz)
       < scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                    prx,pry,prz,psx,psy,psz,
                                    px,py,pz);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
has_smaller_signed_dist_to_planeC3(
     const FT &ppx, const FT &ppy, const FT &ppz,
     const FT &pqx, const FT &pqy, const FT &pqz,
     const FT &prx, const FT &pry, const FT &prz,
     const FT &px, const FT &py, const FT &pz,
     const FT &qx, const FT &qy, const FT &qz)
{
  return scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                    prx,pry,prz,psx,psy,psz,
                                    px,py,pz)
       < scaled_distance_to_planeC3(ppx,ppy,ppz,pqx,pqy,pqz,
                                    prx,pry,prz,psx,psy,psz,
                                    qx,qy,qz);
}


CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_PREDICATES_ON_FTC3_H
#include <CGAL/Arithmetic_filter/predicates_on_ftC3.h>
#endif // CGAL_ARITHMETIC_FILTER_PREDICATES_ON_FTC3_H
#endif

#endif // CGAL_PREDICATES_KERNEL_FTC3_H
