// Copyright (c) 2000, 2016
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Herve Bronnimann, Sylvain Pion, Olivier Devillers
//                

#ifndef CGAL_PREDICATES_KERNEL_FTC2_H
#define CGAL_PREDICATES_KERNEL_FTC2_H

#include <CGAL/algorithm.h>
#include <CGAL/number_utils.h>
#include <CGAL/predicates/sign_of_determinant.h>
#include <CGAL/constructions/kernel_ftC2.h>

namespace CGAL {

template < class FT >
inline
typename Equal_to<FT>::result_type
parallelC2(const FT &l1a, const FT &l1b,
           const FT &l2a, const FT &l2b)
{
    return sign_of_determinant(l1a, l1b, l2a, l2b) == ZERO;
}

template < class FT >
typename Equal_to<FT>::result_type
parallelC2(const FT &s1sx, const FT &s1sy,
           const FT &s1tx, const FT &s1ty,
           const FT &s2sx, const FT &s2sy,
           const FT &s2tx, const FT &s2ty)
{
    return sign_of_determinant(s1tx - s1sx, s1ty - s1sy,
                                  s2tx - s2sx, s2ty - s2sy) == ZERO;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Equal_to<FT>::result_type
equal_lineC2(const FT &l1a, const FT &l1b, const FT &l1c,
             const FT &l2a, const FT &l2b, const FT &l2c)
{
    if (sign_of_determinant(l1a, l1b, l2a, l2b) != ZERO)
        return false; // Not parallel.
    typename Sgn<FT>::result_type s1a = CGAL_NTS sign(l1a);
    if (s1a != ZERO)
        return s1a == CGAL_NTS sign(l2a)
	    && sign_of_determinant(l1a, l1c, l2a, l2c) == ZERO;
    return CGAL_NTS sign(l1b) == CGAL_NTS sign(l2b)
	&& sign_of_determinant(l1b, l1c, l2b, l2c) == ZERO;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_xC2(const FT &px,
            const FT &la, const FT &lb, const FT &lc,
            const FT &ha, const FT &hb, const FT &hc)
{
  // The abscissa of the intersection point is num/den.
  FT num = determinant( lb, lc, hb, hc);
  FT den = determinant( la, lb, ha, hb);
  typename Sgn<FT>::result_type s = CGAL_NTS sign(den);
  CGAL_kernel_assertion( s != ZERO );
  return s * CGAL_NTS compare(px * den, num);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_xC2(const FT &la, const FT &lb, const FT &lc,
            const FT &h1a, const FT &h1b, const FT &h1c,
            const FT &h2a, const FT &h2b, const FT &h2c)
{
  /*
  FT num1 = determinant( lb, lc, h1b, h1c);
  FT den1 = determinant( la, lb, h1a, h1b);
  FT num2 = determinant( lb, lc, h2b, h2c);
  FT den2 = determinant( la, lb, h2a, h2b);
  Sign s = Sign (CGAL_NTS sign(den1) * CGAL_NTS sign(den2));
  CGAL_kernel_assertion( s != ZERO );
  return s * sign_of_determinant(num1, num2, den1, den2);
  */
  FT num1 = determinant( la, lc, h1a, h1c);
  FT num2 = determinant( la, lc, h2a, h2c);
  FT num  = determinant(h1a,h1c,h2a,h2c)*lb
            + determinant(num1,num2,h1b,h2b);
  FT den1 = determinant( la, lb, h1a, h1b);
  FT den2 = determinant( la, lb, h2a, h2b);
  return CGAL_NTS sign(lb) *
         CGAL_NTS sign(num) *
         CGAL_NTS sign(den1) *
         CGAL_NTS sign(den2);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_xC2(const FT &l1a, const FT &l1b, const FT &l1c,
            const FT &h1a, const FT &h1b, const FT &h1c,
            const FT &l2a, const FT &l2b, const FT &l2c,
            const FT &h2a, const FT &h2b, const FT &h2c)
{
  FT num1 = determinant( l1b, l1c, h1b, h1c);
  FT den1 = determinant( l1a, l1b, h1a, h1b);
  FT num2 = determinant( l2b, l2c, h2b, h2c);
  FT den2 = determinant( l2a, l2b, h2a, h2b);
  typename Sgn<FT>::result_type s = CGAL_NTS sign(den1) * CGAL_NTS sign(den2);
  CGAL_kernel_assertion( s != ZERO );
  return s * sign_of_determinant(num1, num2, den1, den2);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_y_at_xC2(const FT &px, const FT &py,
                 const FT &la, const FT &lb, const FT &lc)
{
  typename Sgn<FT>::result_type s = CGAL_NTS sign(lb);
  CGAL_kernel_assertion( s != ZERO );
  return s * CGAL_NTS sign(la*px + lb*py + lc);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_y_at_xC2(const FT &px,
                 const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c)
{
  typename Sgn<FT>::result_type s = CGAL_NTS sign(l1b) * CGAL_NTS sign(l2b);
  CGAL_kernel_assertion( s != ZERO );
  return s * sign_of_determinant<FT>(l2a*px+l2c, l2b,
                                        l1a*px+l1c, l1b);
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Compare<FT>::result_type
compare_y_at_xC2(const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c,
                 const FT &ha,  const FT &hb,  const FT &hc)
{
  typename Sgn<FT>::result_type s = CGAL_NTS sign(hb) *
                                    sign_of_determinant(l1a, l1b, l2a, l2b);
  CGAL_kernel_assertion( s != ZERO );
  return s * sign_of_determinant(l1a, l1b, l1c,
                                    l2a, l2b, l2c,
                                    ha,  hb,  hc);
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Compare<FT>::result_type
compare_y_at_xC2(const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c,
                 const FT &h1a, const FT &h1b, const FT &h1c,
                 const FT &h2a, const FT &h2b, const FT &h2c)
{
  // The abscissa of the intersection point is num/den.
  FT num = determinant( l1b, l1c, l2b, l2c);
  FT den = determinant( l1a, l1b, l2a, l2b);
  typename Sgn<FT>::result_type s = CGAL_NTS sign(h1b) *
                                    CGAL_NTS sign(h2b) *
                                    CGAL_NTS sign(den);
  CGAL_kernel_assertion( s != ZERO );
  return s * sign_of_determinant<FT>(h2a*num+h2c*den, h2b,
                                        h1a*num+h1c*den, h1b);
}

// forward-declaration of orientationC2, used in compare_y_at_xC2
template < class FT >
inline
typename Same_uncertainty_nt<Orientation, FT>::type
orientationC2(const FT &px, const FT &py,
              const FT &qx, const FT &qy,
              const FT &rx, const FT &ry);

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Compare<FT>::result_type
compare_y_at_xC2(const FT &px, const FT &py,
                 const FT &ssx, const FT &ssy,
                 const FT &stx, const FT &sty)
{
    // compares the y-coordinates of p and the vertical projection of p on s.
    // Precondition : p is in the x-range of s.

    CGAL_kernel_precondition(are_ordered(ssx, px, stx));

    if (ssx < stx)
	return orientationC2(px, py, ssx, ssy, stx, sty);
    else if (ssx > stx)
	return orientationC2(px, py, stx, sty, ssx, ssy);
    else {
	if (py < (CGAL::min)(sty, ssy))
	    return SMALLER;
	if (py > (CGAL::max)(sty, ssy))
	    return LARGER;
	return EQUAL;
    }
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Compare<FT>::result_type
compare_y_at_x_segment_C2(const FT &px,
                          const FT &s1sx, const FT &s1sy,
                          const FT &s1tx, const FT &s1ty,
                          const FT &s2sx, const FT &s2sy,
                          const FT &s2tx, const FT &s2ty)
{
    // compares the y-coordinates of the vertical projections of p on s1 and s2
    // Precondition : p is in the x-range of s1 and s2.
    // - if one or two segments are vertical :
    //   - if the segments intersect, return EQUAL
    //   - if not, return the obvious SMALLER/LARGER.

    CGAL_kernel_precondition(are_ordered(s1sx, px, s1tx));
    CGAL_kernel_precondition(are_ordered(s2sx, px, s2tx));

    if (s1sx != s1tx && s2sx != s2tx) {
	FT s1stx = s1sx-s1tx;
	FT s2stx = s2sx-s2tx;

	return CGAL_NTS compare(s1sx, s1tx) *
	       CGAL_NTS compare(s2sx, s2tx) *
	       CGAL_NTS compare(-(s1sx-px)*(s1sy-s1ty)*s2stx,
		                (s2sy-s1sy)*s2stx*s1stx
		                -(s2sx-px)*(s2sy-s2ty)*s1stx );
    }
    else {
	if (s1sx == s1tx) { // s1 is vertical
	    typename Compare<FT>::result_type c1, c2;
	    c1 = compare_y_at_xC2(px, s1sy, s2sx, s2sy, s2tx, s2ty);
	    c2 = compare_y_at_xC2(px, s1ty, s2sx, s2sy, s2tx, s2ty);
	    if (c1 == c2)
		return c1;
	    return EQUAL;
	}
	// s2 is vertical
	typename Compare<FT>::result_type c3, c4;
	c3 = compare_y_at_xC2(px, s2sy, s1sx, s1sy, s1tx, s1ty);
	c4 = compare_y_at_xC2(px, s2ty, s1sx, s1sy, s1tx, s1ty);
	if (c3 == c4)
	    return -c3;
	return EQUAL;
    }
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Equal_to<FT>::result_type
equal_directionC2(const FT &dx1, const FT &dy1,
                  const FT &dx2, const FT &dy2) 
{
  return CGAL_NTS sign(dx1) == CGAL_NTS sign(dx2)
      && CGAL_NTS sign(dy1) == CGAL_NTS sign(dy2)
      && sign_of_determinant(dx1, dy1, dx2, dy2) == ZERO;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_angle_with_x_axisC2(const FT &dx1, const FT &dy1,
                            const FT &dx2, const FT &dy2) 
{
  // angles are in [-pi,pi], and the angle between Ox and d1 is compared
  // with the angle between Ox and d2
  int quadrant_1 = (dx1 >= 0) ? (dy1 >= 0 ? 1 : 4)
                              : (dy1 >= 0 ? 2 : 3);
  int quadrant_2 = (dx2 >= 0) ? (dy2 >= 0 ? 1 : 4)
                              : (dy2 >= 0 ? 2 : 3);
  // We can't use CGAL_NTS compare(quadrant_1,quadrant_2) because in case
  // of tie, we need additional computation
  if (quadrant_1 > quadrant_2)
    return LARGER;
  else if (quadrant_1 < quadrant_2)
    return SMALLER;
  return -sign_of_determinant(dx1,dy1,dx2,dy2);
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_slopesC2(const FT &l1a, const FT &l1b, const FT &l2a, const FT &l2b) 
{
   typedef typename Compare<FT>::result_type result_type;

   if (CGAL_NTS is_zero(l1a))  // l1 is horizontal
    return CGAL_NTS is_zero(l2b) ? result_type(SMALLER)
	                         : CGAL_NTS sign(l2a) * CGAL_NTS sign(l2b);
   if (CGAL_NTS is_zero(l2a)) // l2 is horizontal
    return CGAL_NTS is_zero(l1b) ? result_type(LARGER)
	                         : - CGAL_NTS sign(l1a) * CGAL_NTS sign(l1b);
   if (CGAL_NTS is_zero(l1b)) return CGAL_NTS is_zero(l2b) ? EQUAL : LARGER;
   if (CGAL_NTS is_zero(l2b)) return SMALLER;
   result_type l1_sign = - CGAL_NTS sign(l1a) * CGAL_NTS sign(l1b);
   result_type l2_sign = - CGAL_NTS sign(l2a) * CGAL_NTS sign(l2b);

   if (l1_sign < l2_sign) return SMALLER;
   if (l1_sign > l2_sign) return LARGER;

   if (l1_sign > ZERO)
     return CGAL_NTS compare ( CGAL_NTS abs(l1a * l2b),
			       CGAL_NTS abs(l2a * l1b) );

   return CGAL_NTS compare ( CGAL_NTS abs(l2a * l1b),
			     CGAL_NTS abs(l1a * l2b) );
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Compare<FT>::result_type
compare_slopesC2(const FT &s1_src_x, const FT &s1_src_y, const FT &s1_tgt_x, 
                 const FT &s1_tgt_y, const FT &s2_src_x, const FT &s2_src_y, 
                 const FT &s2_tgt_x, const FT &s2_tgt_y) 
{
   typedef typename Compare<FT>::result_type  Cmp;
   typedef typename Sgn<FT>::result_type      Sg;
   Cmp cmp_y1 = CGAL_NTS compare(s1_src_y, s1_tgt_y);
   if (cmp_y1 == EQUAL) // horizontal
   {
      Cmp cmp_x2 = CGAL_NTS compare(s2_src_x, s2_tgt_x);

      if (cmp_x2 == EQUAL) return SMALLER;
      return - CGAL_NTS sign(s2_src_y - s2_tgt_y) * CGAL_NTS sign(s2_src_x - s2_tgt_x);
   }

   Cmp cmp_y2 = CGAL_NTS compare(s2_src_y, s2_tgt_y);
   if (cmp_y2 == EQUAL)
   {
      Cmp cmp_x1 = CGAL_NTS compare(s1_src_x, s1_tgt_x);

      if (cmp_x1 == EQUAL) return LARGER;
      return CGAL_NTS sign(s1_src_y - s1_tgt_y) * CGAL_NTS sign(s1_src_x - s1_tgt_x);
   }

   Cmp cmp_x1 = CGAL_NTS compare(s1_src_x, s1_tgt_x);
   Cmp cmp_x2 = CGAL_NTS compare(s2_src_x, s2_tgt_x);

   if (cmp_x1 == EQUAL) return cmp_x2 == EQUAL ? EQUAL : LARGER;

   if (cmp_x2 == EQUAL) return SMALLER;

   FT s1_xdiff = s1_src_x - s1_tgt_x;
   FT s1_ydiff = s1_src_y - s1_tgt_y;
   FT s2_xdiff = s2_src_x - s2_tgt_x;
   FT s2_ydiff = s2_src_y - s2_tgt_y;
   Sg s1_sign = CGAL_NTS sign(s1_ydiff) * CGAL_NTS sign(s1_xdiff);
   Sg s2_sign = CGAL_NTS sign(s2_ydiff) * CGAL_NTS sign(s2_xdiff);

   if (s1_sign < s2_sign) return SMALLER;
   if (s1_sign > s2_sign) return LARGER;

   if (s1_sign > ZERO)
     return CGAL_NTS compare( CGAL_NTS abs(s1_ydiff * s2_xdiff),
                              CGAL_NTS abs(s2_ydiff * s1_xdiff));

   return CGAL_NTS compare( CGAL_NTS abs(s2_ydiff * s1_xdiff),
                            CGAL_NTS abs(s1_ydiff * s2_xdiff));
}


#if 0
// Unused, undocumented, un-functorized.
template < class FT >
inline
typename Compare<FT>::result_type
compare_deltax_deltayC2(const FT &px, const FT &qx,
                        const FT &ry, const FT &sy)
{
  return CGAL_NTS compare(CGAL_NTS abs(px-qx), CGAL_NTS abs(ry-sy));
}
#endif

template < class FT >
inline
typename Compare<FT>::result_type
compare_lexicographically_xyC2(const FT &px, const FT &py,
                               const FT &qx, const FT &qy)
{
  typename Compare<FT>::result_type c = CGAL_NTS compare(px,qx);
  return (c != EQUAL) ? c : CGAL_NTS compare(py,qy);
}

template < class FT >
inline
typename Same_uncertainty_nt<Orientation, FT>::type
orientationC2(const FT &px, const FT &py,
              const FT &qx, const FT &qy,
              const FT &rx, const FT &ry)
{
  return sign_of_determinant(qx-px, qy-py, rx-px, ry-py);
}

template < class FT >
inline
typename Same_uncertainty_nt<Orientation, FT>::type
orientationC2(const FT &ux, const FT &uy, const FT &vx, const FT &vy)
{
  return sign_of_determinant(ux, uy, vx, vy);
}

template < class FT >
inline
typename Same_uncertainty_nt<Angle, FT>::type
angleC2(const FT &ux, const FT &uy,
        const FT &vx, const FT &vy)
{
  return enum_cast<Angle>(CGAL_NTS sign(ux*vx + uy*vy));
}

template < class FT >
inline
typename Same_uncertainty_nt<Angle, FT>::type
angleC2(const FT &px, const FT &py,
        const FT &qx, const FT &qy,
        const FT &rx, const FT &ry)
{
  return enum_cast<Angle>(CGAL_NTS sign((px-qx)*(rx-qx)+(py-qy)*(ry-qy)));
}

template < class FT >
inline
typename Same_uncertainty_nt<Angle, FT>::type
angleC2(const FT &px, const FT &py,
        const FT &qx, const FT &qy,
        const FT &rx, const FT &ry,
        const FT &sx, const FT &sy)
{
  return enum_cast<Angle>(CGAL_NTS sign((px-qx)*(rx-sx)+(py-qy)*(ry-sy)));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Equal_to<FT>::result_type
collinear_are_ordered_along_lineC2(const FT &px, const FT &py,
                                   const FT &qx, const FT &qy,
                                   const FT &rx, const FT &ry)
{
  if (px < qx) return !(rx < qx);
  if (qx < px) return !(qx < rx);
  if (py < qy) return !(ry < qy);
  if (qy < py) return !(qy < ry);
  return true; // p==q
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
typename Equal_to<FT>::result_type
collinear_are_strictly_ordered_along_lineC2(const FT &px, const FT &py,
                                            const FT &qx, const FT &qy,
                                            const FT &rx, const FT &ry)
{
  if (px < qx) return (qx < rx);
  if (qx < px) return (rx < qx);
  if (py < qy) return (qy < ry);
  if (qy < py) return (ry < qy);
  return false;
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Same_uncertainty_nt<Oriented_side, FT>::type
side_of_oriented_circleC2(const FT &px, const FT &py,
                          const FT &qx, const FT &qy,
                          const FT &rx, const FT &ry,
                          const FT &tx, const FT &ty)
{
  //  sign_of_determinant(px, py, px*px + py*py, 1,
  //                         qx, qy, qx*qx + qy*qy, 1,
  //                         rx, ry, rx*rx + ry*ry, 1,
  //                         tx, ty, tx*tx + ty*ty, 1);
  // We first translate so that p is the new origin.
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT tpx = tx-px;
  FT tpy = ty-py;
// The usual 3x3 formula can be simplified a little bit to a 2x2.
//         - sign_of_determinant(qpx, qpy, square(qpx) + square(qpy),
//                                  rpx, rpy, square(rpx) + square(rpy),
//                                  tpx, tpy, square(tpx) + square(tpy)));
  return sign_of_determinant<FT>( qpx*tpy - qpy*tpx, tpx*(tx-qx) + tpy*(ty-qy),
                                  qpx*rpy - qpy*rpx, rpx*(rx-qx) + rpy*(ry-qy));
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Same_uncertainty_nt<Bounded_side, FT>::type
side_of_bounded_circleC2(const FT &px, const FT &py,
                         const FT &qx, const FT &qy,
                         const FT &rx, const FT &ry,
                         const FT &tx, const FT &ty)
{
  return enum_cast<Bounded_side>( side_of_oriented_circleC2(px,py,qx,qy,rx,ry,tx,ty)
                                 * orientationC2(px,py,qx,qy,rx,ry) );
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
typename Same_uncertainty_nt<Bounded_side, FT>::type
side_of_bounded_circleC2(const FT &px, const FT &py,
                         const FT &qx, const FT &qy,
                         const FT &tx, const FT &ty)
{
  // Returns whether T lies inside or outside the circle which diameter is PQ.
  return enum_cast<Bounded_side>(
                      CGAL_NTS compare((tx-px)*(qx-tx), (ty-py)*(ty-qy)) );
}

template < class FT >
inline
typename Compare<FT>::result_type
cmp_dist_to_pointC2(const FT &px, const FT &py,
                    const FT &qx, const FT &qy,
                    const FT &rx, const FT &ry)
{
  return CGAL_NTS compare(squared_distanceC2(px,py,qx,qy),
                          squared_distanceC2(px,py,rx,ry));
}

template < class FT >
inline
typename Equal_to<FT>::result_type
has_larger_dist_to_pointC2(const FT &px, const FT &py,
                           const FT &qx, const FT &qy,
                           const FT &rx, const FT &ry)
{
  return cmp_dist_to_pointC2(px,py,qx,qy,rx,ry) == LARGER;
}

template < class FT >
inline
typename Equal_to<FT>::result_type
has_smaller_dist_to_pointC2(const FT &px, const FT &py,
                            const FT &qx, const FT &qy,
                            const FT &rx, const FT &ry)
{
  return cmp_dist_to_pointC2(px,py,qx,qy,rx,ry) == SMALLER;
}

template < class FT >
inline
typename Compare<FT>::result_type
cmp_signed_dist_to_directionC2(const FT &la, const FT &lb,
                               const FT &px, const FT &py,
                               const FT &qx, const FT &qy)
{
  return CGAL_NTS compare(scaled_distance_to_directionC2(la,lb,px,py),
                          scaled_distance_to_directionC2(la,lb,qx,qy));
}

template < class FT >
inline
typename Equal_to<FT>::result_type
has_larger_signed_dist_to_directionC2(const FT &la, const FT &lb,
                                      const FT &px, const FT &py,
                                      const FT &qx, const FT &qy)
{
  return cmp_signed_dist_to_directionC2(la,lb,px,py,qx,qy) == LARGER;
}

template < class FT >
inline
typename Equal_to<FT>::result_type
has_smaller_signed_dist_to_directionC2(const FT &la, const FT &lb,
                                       const FT &px, const FT &py,
                                       const FT &qx, const FT &qy)
{
  return cmp_signed_dist_to_directionC2(la,lb,px,py,qx,qy) == SMALLER;
}

template <class FT>
inline
typename Compare<FT>::result_type
cmp_signed_dist_to_lineC2(const FT &px, const FT &py,
                          const FT &qx, const FT &qy,
                          const FT &rx, const FT &ry,
                          const FT &sx, const FT &sy)
{
  return CGAL_NTS compare(scaled_distance_to_lineC2(px,py,qx,qy,rx,ry),
                          scaled_distance_to_lineC2(px,py,qx,qy,sx,sy));
}

template <class FT>
inline
typename Equal_to<FT>::result_type
has_larger_signed_dist_to_lineC2(const FT &px, const FT &py,
                                 const FT &qx, const FT &qy,
                                 const FT &rx, const FT &ry,
                                 const FT &sx, const FT &sy)
{
  return cmp_signed_dist_to_lineC2(px,py,qx,qy,rx,ry,sx,sy) == LARGER;
}

template <class FT>
inline
typename Equal_to<FT>::result_type
has_smaller_signed_dist_to_lineC2(const FT &px, const FT &py,
                                  const FT &qx, const FT &qy,
                                  const FT &rx, const FT &ry,
                                  const FT &sx, const FT &sy)
{
  return cmp_signed_dist_to_lineC2(px,py,qx,qy,rx,ry,sx,sy) == SMALLER;
}

template <class FT>
inline
typename Same_uncertainty_nt<Oriented_side, FT>::type
side_of_oriented_lineC2(const FT &a, const FT &b, const FT &c,
                        const FT &x, const FT &y)
{
  return CGAL_NTS sign(a*x+b*y+c);
}

template <class FT>
Comparison_result
compare_power_distanceC2(const FT& px, const FT& py, const FT& pwt,
                         const FT& qx, const FT& qy, const FT& qwt,
                         const FT& rx, const FT& ry)
{
  // returns SMALLER if r is closer to p w.r.t. the power metric
  FT d1 = CGAL_NTS square(rx - px) + CGAL_NTS square(ry - py) - pwt;
  FT d2 = CGAL_NTS square(rx - qx) + CGAL_NTS square(ry - qy) - qwt;
  return CGAL_NTS compare(d1, d2);
}

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
power_side_of_bounded_power_circleC2(const FT &px, const FT &py, const FT &pw,
                                     const FT &qx, const FT &qy, const FT &qw,
                                     const FT &tx, const FT &ty, const FT &tw)
{
  FT dpx = px - qx;
  FT dpy = py - qy;
  FT dtx = tx - qx;
  FT dty = ty - qy;
  FT dpz = CGAL_NTS square(dpx) + CGAL_NTS square(dpy);

  return enum_cast<Bounded_side>
    (CGAL_NTS sign(-(CGAL_NTS square(dtx) + CGAL_NTS square(dty)-tw+qw)*dpz
                   +(dpz-pw+qw)*(dpx*dtx+dpy*dty)));
}

template <class FT>
Oriented_side
power_side_of_oriented_power_circleC2(const FT &px, const FT &py, const FT &pwt,
                                      const FT &qx, const FT &qy, const FT &qwt,
                                      const FT &rx, const FT &ry, const FT &rwt,
                                      const FT &tx, const FT &ty, const FT &twt)
{
  // Note: maybe this can be further optimized like the usual in_circle() ?

  // We translate the 4 points so that T becomes the origin.
  FT dpx = px - tx;
  FT dpy = py - ty;
  FT dpz = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) - pwt + twt;
  FT dqx = qx - tx;
  FT dqy = qy - ty;
  FT dqz = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - qwt + twt;
  FT drx = rx - tx;
  FT dry = ry - ty;
  FT drz = CGAL_NTS square(drx) + CGAL_NTS square(dry) - rwt + twt;

  return sign_of_determinant(dpx, dpy, dpz,
                             dqx, dqy, dqz,
                             drx, dry, drz);
}

template <class FT>
Oriented_side
power_side_of_oriented_power_circleC2(const FT &px, const FT &py, const FT &pwt,
                                      const FT &qx, const FT &qy, const FT &qwt,
                                      const FT &tx, const FT &ty, const FT &twt)
{
  // Same translation as above.
  FT dpx = px - tx;
  FT dpy = py - ty;
  FT dpz = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) - pwt + twt;
  FT dqx = qx - tx;
  FT dqy = qy - ty;
  FT dqz = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - qwt + twt;

  // We do an orthogonal projection on the (x) axis, if possible.
  Comparison_result cmpx = CGAL_NTS compare(px, qx);
  if (cmpx != EQUAL)
    return cmpx * sign_of_determinant(dpx, dpz, dqx, dqz);

  // If not possible, then on the (y) axis.
  Comparison_result cmpy = CGAL_NTS compare(py, qy);
  return cmpy * sign_of_determinant(dpy, dpz, dqy, dqz);
}

} //namespace CGAL

#endif  // CGAL_PREDICATES_KERNEL_FTC2_H
