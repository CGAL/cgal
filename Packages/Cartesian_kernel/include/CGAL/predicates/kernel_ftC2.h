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
// file          : include/CGAL/predicates/kernel_ftC2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_PREDICATES_KERNEL_FTC2_H
#define CGAL_PREDICATES_KERNEL_FTC2_H

#include <CGAL/number_utils.h>
#include <CGAL/predicates/sign_of_determinant.h>
#include <CGAL/constructions/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
equal_lineC2(const FT &l1a, const FT &l1b, const FT &l1c,
             const FT &l2a, const FT &l2b, const FT &l2c)
{
    if (sign_of_determinant2x2(l1a, l1b, l2a, l2b) != ZERO)
        return false; // Not parallel.
    CGAL::Sign s1a = CGAL_NTS sign(l1a);
    if (s1a != ZERO)
        return s1a == CGAL_NTS sign(l2a)
	    && sign_of_determinant2x2(l1a, l1c, l2a, l2c) == ZERO;
    return CGAL_NTS sign(l1b) == CGAL_NTS sign(l2b)
	&& sign_of_determinant2x2(l1b, l1c, l2b, l2c) == ZERO;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_xC2(const FT &px,
            const FT &la, const FT &lb, const FT &lc,
            const FT &ha, const FT &hb, const FT &hc)
{
  // The abscissa of the intersection point is num/den.
  FT num = det2x2_by_formula( lb, lc, hb, hc);
  FT den = det2x2_by_formula( la, lb, ha, hb);
  Sign s = CGAL_NTS sign(den);
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result( s * CGAL_NTS compare( px * den, num) );
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_xC2(const FT &la, const FT &lb, const FT &lc,
            const FT &h1a, const FT &h1b, const FT &h1c,
            const FT &h2a, const FT &h2b, const FT &h2c)
{
  /*
  FT num1 = det2x2_by_formula( lb, lc, h1b, h1c);
  FT den1 = det2x2_by_formula( la, lb, h1a, h1b);
  FT num2 = det2x2_by_formula( lb, lc, h2b, h2c);
  FT den2 = det2x2_by_formula( la, lb, h2a, h2b);
  Sign s = Sign (CGAL_NTS sign(den1) * CGAL_NTS sign(den2));
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result( s * sign_of_determinant2x2(num1, 
                                                       num2, den1, den2));
  */
  FT num1 = det2x2_by_formula( la, lc, h1a, h1c);
  FT num2 = det2x2_by_formula( la, lc, h2a, h2c);
  FT num  = det2x2_by_formula(h1a,h1c,h2a,h2c)*lb
            + det2x2_by_formula(num1,num2,h1b,h2b);
  FT den1 = det2x2_by_formula( la, lb, h1a, h1b);
  FT den2 = det2x2_by_formula( la, lb, h2a, h2b);
  return Comparison_result( CGAL_NTS sign(lb) * CGAL_NTS sign(num) *
                            CGAL_NTS sign(den1) * CGAL_NTS sign(den2));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_xC2(const FT &l1a, const FT &l1b, const FT &l1c,
            const FT &h1a, const FT &h1b, const FT &h1c,
            const FT &l2a, const FT &l2b, const FT &l2c,
            const FT &h2a, const FT &h2b, const FT &h2c)
{
  FT num1 = det2x2_by_formula( l1b, l1c, h1b, h1c);
  FT den1 = det2x2_by_formula( l1a, l1b, h1a, h1b);
  FT num2 = det2x2_by_formula( l2b, l2c, h2b, h2c);
  FT den2 = det2x2_by_formula( l2a, l2b, h2a, h2b);
  Sign s = Sign (CGAL_NTS sign(den1) * CGAL_NTS sign(den2));
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result( s * sign_of_determinant2x2(num1, num2,
			                               den1, den2));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_y_at_xC2(const FT &px, const FT &py,
                 const FT &la, const FT &lb, const FT &lc)
{
  Sign s = CGAL_NTS sign(lb);
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result (s * CGAL_NTS sign(la*px + lb*py + lc));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_y_at_xC2(const FT &px,
                 const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c)
{
  Sign s = Sign (CGAL_NTS sign(l1b) * CGAL_NTS sign(l2b));
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result ( s * sign_of_determinant2x2(l2a*px+l2c, l2b,
                                                        l1a*px+l1c, l1b));
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
compare_y_at_xC2(const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c,
                 const FT &ha,  const FT &hb,  const FT &hc)
{
  Sign s = Sign (sign_of_determinant2x2(l1a, l1b, l2a, l2b) *
		 CGAL_NTS sign(hb));
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result( s * sign_of_determinant3x3(l1a, l1b, l1c,
                                                       l2a, l2b, l2c,
                                                       ha,  hb,  hc));
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
compare_y_at_xC2(const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c,
                 const FT &h1a, const FT &h1b, const FT &h1c,
                 const FT &h2a, const FT &h2b, const FT &h2c)
{
  // The abscissa of the intersection point is num/den.
  FT num = det2x2_by_formula( l1b, l1c, l2b, l2c);
  FT den = det2x2_by_formula( l1a, l1b, l2a, l2b);
  Sign s = Sign (CGAL_NTS sign(h1b) * CGAL_NTS sign(h2b) * CGAL_NTS sign(den));
  CGAL_kernel_assertion( s != ZERO );
  return Comparison_result ( s * sign_of_determinant2x2(h2a*num+h2c*den, h2b,
                                                        h1a*num+h1c*den, h1b));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
equal_directionC2(const FT &dx1, const FT &dy1,
                  const FT &dx2, const FT &dy2) 
{
  return CGAL_NTS sign(dx1) == CGAL_NTS sign(dx2)
      && CGAL_NTS sign(dy1) == CGAL_NTS sign(dy2)
      && sign_of_determinant2x2(dx1, dy1, dx2, dy2) == ZERO;
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_angle_with_x_axisC2(const FT &dx1, const FT &dy1,
                            const FT &dx2, const FT &dy2) 
{
  // angles are in [-pi,pi], and the angle between Ox and d1 is compared
  // with the angle between Ox and d2
  int quadrant_1 = (dx1 >= FT(0)) ? ((dy1 >= FT(0))?1:4)
                                  : ((dy1 >= FT(0))?2:3);
  int quadrant_2 = (dx2 >= FT(0)) ? ((dy2 >= FT(0))?1:4)
                                  : ((dy2 >= FT(0))?2:3);
  // We can't use CGAL_NTS compare(quadrant_1,quadrant_2) because in case
  // of tie, we need additional computation
  if (quadrant_1 > quadrant_2)
    return LARGER;
  else if (quadrant_1 < quadrant_2)
    return SMALLER;
  return Comparison_result(-sign_of_determinant2x2(dx1,dy1,dx2,dy2));
}

template < class FT >
inline
Comparison_result
compare_deltax_deltayC2(const FT &px, const FT &qx,
                        const FT &ry, const FT &sy)
{
  return CGAL_NTS compare(CGAL_NTS abs(px-qx), CGAL_NTS abs(ry-sy));
}

template < class FT >
/*CGAL_NO_FILTER*/
inline
Comparison_result
compare_lexicographically_xyC2(const FT &px, const FT &py,
                               const FT &qx, const FT &qy)
{
  Comparison_result c = CGAL_NTS compare(px,qx);
  return (c != EQUAL) ? c : CGAL_NTS compare(py,qy);
}

template < class FT >
inline
Orientation
orientationC2(const FT &px, const FT &py,
              const FT &qx, const FT &qy,
              const FT &rx, const FT &ry)
{
  return Orientation (sign_of_determinant2x2(qx-px, qy-py,
                                             rx-px, ry-py));
}

template < class FT >
inline
Angle
angleC2(const FT &px, const FT &py,
        const FT &qx, const FT &qy,
        const FT &rx, const FT &ry)
{
  return (Angle) CGAL_NTS sign ((px-qx)*(rx-qx)+(py-qy)*(ry-qy));
}

template < class FT >
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
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
/*CGAL_NO_FILTER*/
CGAL_KERNEL_MEDIUM_INLINE
bool
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
Oriented_side
side_of_oriented_circleC2(const FT &px, const FT &py,
                          const FT &qx, const FT &qy,
                          const FT &rx, const FT &ry,
                          const FT &tx, const FT &ty)
{
  // Oriented_side(
  //  sign_of_determinant4x4(px, py, px*px + py*py, 1,
  //                         qx, qy, qx*qx + qy*qy, 1,
  //                         rx, ry, rx*rx + ry*ry, 1,
  //                         tx, ty, tx*tx + ty*ty, 1));
  // We first translate so that p is the new origin.
  FT qpx = qx-px;
  FT qpy = qy-py;
  FT rpx = rx-px;
  FT rpy = ry-py;
  FT tpx = tx-px;
  FT tpy = ty-py;
// The usual 3x3 formula can be simplified a little bit to a 2x2.
//         - sign_of_determinant3x3(qpx, qpy, square(qpx) + square(qpy),
//                                  rpx, rpy, square(rpx) + square(rpy),
//                                  tpx, tpy, square(tpx) + square(tpy)));
  return Oriented_side(sign_of_determinant2x2(
                             qpx*tpy - qpy*tpx, tpx*(tx-qx) + tpy*(ty-qy),
                             qpx*rpy - qpy*rpx, rpx*(rx-qx) + rpy*(ry-qy)));
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
side_of_bounded_circleC2(const FT &px, const FT &py,
                         const FT &qx, const FT &qy,
                         const FT &rx, const FT &ry,
                         const FT &tx, const FT &ty)
{
  return Bounded_side( side_of_oriented_circleC2(px,py,qx,qy,rx,ry,tx,ty)
                                 * orientationC2(px,py,qx,qy,rx,ry) );
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
side_of_bounded_circleC2(const FT &px, const FT &py,
                         const FT &qx, const FT &qy,
                         const FT &tx, const FT &ty)
{
  // Returns whether T lies inside or outside the circle which diameter is PQ.
  return Bounded_side( CGAL_NTS compare((tx-px)*(qx-tx), (ty-py)*(ty-qy)) );
}

template < class FT >
inline
Comparison_result
cmp_dist_to_pointC2(const FT &px, const FT &py,
                    const FT &qx, const FT &qy,
                    const FT &rx, const FT &ry)
{
  return CGAL_NTS compare(squared_distanceC2(px,py,qx,qy),
                          squared_distanceC2(px,py,rx,ry));
}

template < class FT >
/*CGAL_NO_FILTER*/
inline
bool
has_larger_dist_to_pointC2(const FT &px, const FT &py,
                           const FT &qx, const FT &qy,
                           const FT &rx, const FT &ry)
{
  return cmp_dist_to_pointC2(px,py,qx,qy,rx,ry) == LARGER;
}

template < class FT >
/*CGAL_NO_FILTER*/
inline
bool
has_smaller_dist_to_pointC2(const FT &px, const FT &py,
                            const FT &qx, const FT &qy,
                            const FT &rx, const FT &ry)
{
  return cmp_dist_to_pointC2(px,py,qx,qy,rx,ry) == SMALLER;
}

template < class FT >
inline
Comparison_result
cmp_signed_dist_to_directionC2(const FT &la, const FT &lb,
                               const FT &px, const FT &py,
                               const FT &qx, const FT &qy)
{
  return CGAL_NTS compare(scaled_distance_to_directionC2(la,lb,px,py),
                          scaled_distance_to_directionC2(la,lb,qx,qy));
}

template < class FT >
/*CGAL_NO_FILTER*/
inline
bool
has_larger_signed_dist_to_directionC2(const FT &la, const FT &lb,
                                      const FT &px, const FT &py,
                                      const FT &qx, const FT &qy)
{
  return cmp_signed_dist_to_directionC2(la,lb,px,py,qx,qy) == LARGER;
}

template < class FT >
/*CGAL_NO_FILTER*/
inline
bool
has_smaller_signed_dist_to_directionC2(const FT &la, const FT &lb,
                                       const FT &px, const FT &py,
                                       const FT &qx, const FT &qy)
{
  return cmp_signed_dist_to_directionC2(la,lb,px,py,qx,qy) == SMALLER;
}

template <class FT>
inline
Comparison_result
cmp_signed_dist_to_lineC2(const FT &px, const FT &py,
                          const FT &qx, const FT &qy,
                          const FT &rx, const FT &ry,
                          const FT &sx, const FT &sy)
{
  return CGAL_NTS compare(scaled_distance_to_lineC2(px,py,qx,qy,rx,ry),
                          scaled_distance_to_lineC2(px,py,qx,qy,sx,sy));
}

template <class FT>
/*CGAL_NO_FILTER*/
inline
bool
has_larger_signed_dist_to_lineC2(const FT &px, const FT &py,
                                 const FT &qx, const FT &qy,
                                 const FT &rx, const FT &ry,
                                 const FT &sx, const FT &sy)
{
  return cmp_signed_dist_to_lineC2(px,py,qx,qy,rx,ry,sx,sy) == LARGER;
}

template <class FT>
/*CGAL_NO_FILTER*/
inline
bool
has_smaller_signed_dist_to_lineC2(const FT &px, const FT &py,
                                  const FT &qx, const FT &qy,
                                  const FT &rx, const FT &ry,
                                  const FT &sx, const FT &sy)
{
  return cmp_signed_dist_to_lineC2(px,py,qx,qy,rx,ry,sx,sy) == SMALLER;
}

template <class FT>
inline
Oriented_side
side_of_oriented_lineC2(const FT &a, const FT &b, const FT &c,
                        const FT &x, const FT &y)
{
  return Oriented_side(CGAL_NTS sign(a*x+b*y+c));
}

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/predicates/kernel_ftC2.h>
#endif

#endif  // CGAL_PREDICATES_KERNEL_FTC2_H
