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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr, Herve Bronnimann, Sylvain Pion

#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC2_H
#define CGAL_CONSTRUCTIONS_KERNEL_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>

namespace CGAL {

template < class FT >
CGAL_KERNEL_INLINE
void
midpointC2( const FT &px, const FT &py,
            const FT &qx, const FT &qy,
            FT &x, FT &y )
{
  x = (px+qx) / 2;
  y = (py+qy) / 2;
}

template < class FT >
CGAL_KERNEL_LARGE_INLINE
void
circumcenter_translateC2(const FT &dqx, const FT &dqy,
                         const FT &drx, const FT &dry,
                               FT &dcx,       FT &dcy)
{
  // Given 3 points P, Q, R, this function takes as input:
  // qx-px, qy-py, rx-px, ry-py.  And returns cx-px, cy-py,
  // where (cx, cy) are the coordinates of the circumcenter C.

  // What we do is intersect the bisectors.
  FT r2 = CGAL_NTS square(drx) + CGAL_NTS square(dry);
  FT q2 = CGAL_NTS square(dqx) + CGAL_NTS square(dqy);
  FT den = 2 * determinant(dqx, dqy, drx, dry);

  // The 3 points aren't collinear.
  // Hopefully, this is already checked at the upper level.
  CGAL_kernel_assertion ( ! CGAL_NTS is_zero(den) );

  // One possible optimization here is to precompute 1/den, to avoid one
  // division.  However, we loose precision, and it's maybe not worth it (?).
  dcx =   determinant (dry, dqy, r2, q2) / den;
  dcy = - determinant (drx, dqx, r2, q2) / den;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
void
circumcenterC2( const FT &px, const FT &py,
                const FT &qx, const FT &qy,
                const FT &rx, const FT &ry,
                FT &x, FT &y )
{
  circumcenter_translateC2<FT>(qx-px, qy-py, rx-px, ry-py, x, y);
  x += px;
  y += py;
}

template < class FT >
void
barycenterC2(const FT &p1x, const FT &p1y, const FT &w1,
             const FT &p2x, const FT &p2y,
             FT &x, FT &y)
{
   FT w2 = 1 - w1;
   x = w1 * p1x + w2 * p2x;
   y = w1 * p1y + w2 * p2y;
}

template < class FT >
void
barycenterC2(const FT &p1x, const FT &p1y, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &w2,
             FT &x, FT &y)
{
   FT sum = w1 + w2;
   CGAL_kernel_assertion(sum != 0);
   x = (w1 * p1x + w2 * p2x) / sum;
   y = (w1 * p1y + w2 * p2y) / sum;
}

template < class FT >
void
barycenterC2(const FT &p1x, const FT &p1y, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &w2,
             const FT &p3x, const FT &p3y,
             FT &x, FT &y)
{
   FT w3 = 1 - w1 - w2;
   x = w1 * p1x + w2 * p2x + w3 * p3x;
   y = w1 * p1y + w2 * p2y + w3 * p3y;
}

template < class FT >
void
barycenterC2(const FT &p1x, const FT &p1y, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &w3,
             FT &x, FT &y)
{
   FT sum = w1 + w2 + w3;
   CGAL_kernel_assertion(sum != 0);
   x = (w1 * p1x + w2 * p2x + w3 * p3x) / sum;
   y = (w1 * p1y + w2 * p2y + w3 * p3y) / sum;
}

template < class FT >
void
barycenterC2(const FT &p1x, const FT &p1y, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &w3,
             const FT &p4x, const FT &p4y,
             FT &x, FT &y)
{
   FT w4 = 1 - w1 - w2 - w3;
   x = w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x;
   y = w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y;
}

template < class FT >
void
barycenterC2(const FT &p1x, const FT &p1y, const FT &w1,
             const FT &p2x, const FT &p2y, const FT &w2,
             const FT &p3x, const FT &p3y, const FT &w3,
             const FT &p4x, const FT &p4y, const FT &w4,
             FT &x, FT &y)
{
   FT sum = w1 + w2 + w3 + w4;
   CGAL_kernel_assertion(sum != 0);
   x = (w1 * p1x + w2 * p2x + w3 * p3x + w4 * p4x) / sum;
   y = (w1 * p1y + w2 * p2y + w3 * p3y + w4 * p4y) / sum;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
void
centroidC2( const FT &px, const FT &py,
            const FT &qx, const FT &qy,
            const FT &rx, const FT &ry,
            FT &x, FT &y)
{
   x = (px + qx + rx) / 3;
   y = (py + qy + ry) / 3;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
void
centroidC2( const FT &px, const FT &py,
            const FT &qx, const FT &qy,
            const FT &rx, const FT &ry,
            const FT &sx, const FT &sy,
            FT &x, FT &y)
{
   x = (px + qx + rx + sx) / 4;
   y = (py + qy + ry + sy) / 4;
}

template < class FT >
inline
void
line_from_pointsC2(const FT &px, const FT &py,
                   const FT &qx, const FT &qy,
                   FT &a, FT &b, FT &c) 
{
  // The horizontal and vertical line get a special treatment
  // in order to make the intersection code robust for doubles 
  if(py == qy){
    a = 0 ;
    if(qx > px){
      b = 1;
      c = -py;
    } else if(qx == px){
      b = 0;
      c = 0;
    }else{
      b = -1;
      c = py;
    }
  } else if(qx == px){
    b = 0;
    if(qy > py){
      a = -1;
      c = px;
    } else if (qy == py){
      a = 0;
      c = 0;
    } else {
      a = 1;
      c = -px;
    }
  } else {
    a = py - qy;
    b = qx - px;
    c = -px*a - py*b;
  }
}

template < class FT >
inline
void
line_from_point_directionC2(const FT &px, const FT &py,
                            const FT &dx, const FT &dy,
                            FT &a, FT &b, FT &c) 
{
  a = - dy;
  b = dx;
  c = px*dy - py*dx;
}

template < class FT >
CGAL_KERNEL_INLINE
void
bisector_of_pointsC2(const FT &px, const FT &py,
		     const FT &qx, const FT &qy,
		     FT &a, FT &b, FT& c )
{
  a = 2 * (px - qx);
  b = 2 * (py - qy);
  c = CGAL_NTS square(qx) + CGAL_NTS square(qy) -
      CGAL_NTS square(px) - CGAL_NTS square(py);
}

template < class FT >
CGAL_KERNEL_INLINE
void
bisector_of_linesC2(const FT &pa, const FT &pb, const FT &pc,
		    const FT &qa, const FT &qb, const FT &qc,
		    FT &a, FT &b, FT &c)
{
  // We normalize the equations of the 2 lines, and we then add them.
  FT n1 = CGAL_NTS sqrt(CGAL_NTS square(pa) + CGAL_NTS square(pb));
  FT n2 = CGAL_NTS sqrt(CGAL_NTS square(qa) + CGAL_NTS square(qb));
  a = n2 * pa + n1 * qa;
  b = n2 * pb + n1 * qb;
  c = n2 * pc + n1 * qc;

  // Care must be taken for the case when this produces a degenerate line.
  if (a == 0 && b == 0) {
    a = n2 * pa - n1 * qa;
    b = n2 * pb - n1 * qb;
    c = n2 * pc - n1 * qc;
  }
}

template < class FT >
inline
FT
line_y_at_xC2(const FT &a, const FT &b, const FT &c, const FT &x)
{
  return (-a*x-c) / b;
}

template < class FT > 
inline
void
line_get_pointC2(const FT &a, const FT &b, const FT &c, int i,
                 FT &x, FT &y)
{
  if (CGAL_NTS is_zero(b))
    {
      x = (-b-c)/a + i * b;
      y = 1 - i * a;
    }
  else
    {
      x = 1 + i * b;
      y = -(a+c)/b - i * a;
    }
}

template < class FT > 
inline
void
perpendicular_through_pointC2(const FT &la, const FT &lb,
		              const FT &px, const FT &py,
			      FT &a, FT &b, FT &c)
{
  a = -lb;
  b = la;
  c = lb * px - la * py;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
void
line_project_pointC2(const FT &la, const FT &lb, const FT &lc,
		     const FT &px, const FT &py,
		     FT &x, FT &y)
{
#if 1 // FIXME
  // Original old version
  if (CGAL_NTS is_zero(la)) // horizontal line
  {
    x = px;
    y = -lc/lb;
  }
  else if (CGAL_NTS is_zero(lb)) // vertical line
  {
    x = -lc/la;
    y = py;
  }
  else
  {
    FT ab = la/lb, ba = lb/la, ca = lc/la;
    y = ( -px + ab*py - ca ) / ( ba + ab );
    x = -ba * y - ca;
  }
#else
  // New version, with more multiplications, but less divisions and tests.
  // Let's compare the results of the 2, benchmark them, as well as check
  // the precision with the intervals.
  FT a2 = CGAL_NTS square(la);
  FT b2 = CGAL_NTS square(lb);
  FT d = a2 + b2;
  x = (la * (lb * py - lc) - px * b2) / d;
  y = (lb * (lc - la * px) + py * a2) / d;
#endif
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
squared_radiusC2(const FT &px, const FT &py,
                 const FT &qx, const FT &qy,
                 const FT &rx, const FT &ry,
                 FT &x, FT &y )
{
  circumcenter_translateC2(qx-px, qy-py, rx-px, ry-py, x, y);
  FT r2 = CGAL_NTS square(x) + CGAL_NTS square(y);
  x += px;
  y += py;
  return r2;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
squared_radiusC2(const FT &px, const FT &py,
                 const FT &qx, const FT &qy,
                 const FT &rx, const FT &ry)
{
  FT x, y;
  circumcenter_translateC2<FT>(qx-px, qy-py, rx-px, ry-py, x, y);
  return CGAL_NTS square(x) + CGAL_NTS square(y);
}

template < class FT >
inline
FT
squared_distanceC2( const FT &px, const FT &py,
                    const FT &qx, const FT &qy)
{
  return CGAL_NTS square(px-qx) + CGAL_NTS square(py-qy);
}

template < class FT >
inline
FT
squared_radiusC2(const FT &px, const FT &py,
                 const FT &qx, const FT &qy)
{
  return squared_distanceC2(px, py,qx, qy) / 4;
}

template < class FT >
CGAL_KERNEL_INLINE
FT
scaled_distance_to_lineC2( const FT &la, const FT &lb, const FT &lc,
                           const FT &px, const FT &py)
{
  // for comparisons, use distance_to_directionsC2 instead
  // since lc is irrelevant
  return la*px + lb*py + lc;
}

template < class FT >
CGAL_KERNEL_INLINE
FT
scaled_distance_to_directionC2( const FT &la, const FT &lb,
                                const FT &px, const FT &py)
{
  // scalar product with direction
  return la*px + lb*py;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
scaled_distance_to_lineC2( const FT &px, const FT &py,
                           const FT &qx, const FT &qy,
                           const FT &rx, const FT &ry)
{
  return determinant<FT>(px-rx, py-ry, qx-rx, qy-ry);
}

template < class RT >
void
weighted_circumcenter_translateC2(const RT &dqx, const RT &dqy, const RT &dqw,
                                  const RT &drx, const RT &dry, const RT &drw,
                                  RT &dcx, RT &dcy)
{
  // Given 3 points P, Q, R, this function takes as input:
  // qx-px, qy-py,qw-pw,  rx-px, ry-py, rw-pw.  And returns cx-px, cy-py,
  // where (cx, cy) are the coordinates of the circumcenter C.

  // What we do is intersect the radical axis
  RT r2 = CGAL_NTS square(drx) + CGAL_NTS square(dry) - drw;
  RT q2 = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - dqw;

  RT den = RT(2) * determinant(dqx, dqy, drx, dry);

  // The 3 points aren't collinear.
  // Hopefully, this is already checked at the upper level.
  CGAL_assertion ( den != RT(0) );

  // One possible optimization here is to precompute 1/den, to avoid one
  // division.  However, we loose precision, and it's maybe not worth it (?).
  dcx =   determinant (dry, dqy, r2, q2) / den;
  dcy = - determinant (drx, dqx, r2, q2) / den;
}

//template < class RT >
template < class RT, class We>
void
weighted_circumcenterC2( const RT &px, const RT &py, const We &pw,
                         const RT &qx, const RT &qy, const We &qw,
                         const RT &rx, const RT &ry, const We &rw,
                         RT &x, RT &y )
{
  RT dqw = RT(qw-pw);
  RT drw = RT(rw-pw);

  weighted_circumcenter_translateC2<RT>(qx-px, qy-py, dqw,rx-px, ry-py,drw,x, y);
  x += px;
  y += py;
}

template< class FT >
FT
power_productC2(const FT &px, const FT &py, const FT &pw,
                const FT &qx, const FT &qy, const FT &qw)
{
  // computes the power product of two weighted points
  FT qpx = qx - px;
  FT qpy = qy - py;
  FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy);
  return qp2 - pw - qw;
}

template < class RT , class We>
void
radical_axisC2(const RT &px, const RT &py, const We &pw,
               const RT &qx, const RT &qy, const We &qw,
               RT &a, RT &b, RT& c )
{
  a =  RT(2)*(px - qx);
  b =  RT(2)*(py - qy);
  c = - CGAL_NTS square(px) - CGAL_NTS square(py)
      + CGAL_NTS square(qx) + CGAL_NTS square(qy)
      + RT(pw) - RT(qw);
}

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
squared_radius_orthogonal_circleC2(const FT &px, const FT &py, const FT &pw,
                                   const FT &qx, const FT &qy, const FT &qw,
                                   const FT &rx, const FT &ry, const FT &rw)
{
  FT FT4(4);
  FT dpx = px - rx;
  FT dpy = py - ry;
  FT dqx = qx - rx;
  FT dqy = qy - ry;
  FT dpp = CGAL_NTS square(dpx) + CGAL_NTS square(dpy) - pw + rw;
  FT dqq = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - qw + rw;

  FT det0 = determinant(dpx, dpy, dqx, dqy);
  FT det1 = determinant(dpp, dpy, dqq, dqy);
  FT det2 = determinant(dpx, dpp, dqx, dqq);

  return (CGAL_NTS square(det1) + CGAL_NTS square(det2)) /
                                  (FT4 * CGAL_NTS square(det0)) - rw;
}

template< class FT >
CGAL_KERNEL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonal_circleC2(const FT &px, const FT &py, const FT &pw,
                                            const FT &qx, const FT &qy, const FT &qw)
{
  FT FT4(4);
  FT dpz = CGAL_NTS square(px - qx) + CGAL_NTS square(py - qy);
  return (CGAL_NTS square(dpz - pw + qw) / (FT4 * dpz) - qw);
}

} //namespace CGAL

#endif // CGAL_CONSTRUCTIONS_KERNEL_FTC2_H
