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
// file          : include/CGAL/constructions/kernel_ftC2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr, Hervé Brönnimann, Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC2_H
#define CGAL_CONSTRUCTIONS_KERNEL_FTC2_H

#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_INLINE
void
midpointC2( const FT &px, const FT &py,
            const FT &qx, const FT &qy,
            FT &x, FT &y )
{
  x = (px+qx) / FT(2);
  y = (py+qy) / FT(2);
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
  FT den = FT(2) * det2x2_by_formula(dqx, dqy, drx, dry);

  // The 3 points aren't collinear.
  // Hopefully, this is already checked at the upper level.
  CGAL_kernel_assertion ( ! CGAL_NTS is_zero(den) );

  // One possible optimization here is to precompute 1/den, to avoid one
  // division.  However, we loose precision, and it's maybe not worth it (?).
  dcx =   det2x2_by_formula (dry, dqy, r2, q2) / den;
  dcy = - det2x2_by_formula (drx, dqx, r2, q2) / den;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
void
circumcenterC2( const FT &px, const FT &py,
                const FT &qx, const FT &qy,
                const FT &rx, const FT &ry,
                FT &x, FT &y )
{
  circumcenter_translateC2(qx-px, qy-py, rx-px, ry-py, x, y);
  x += px;
  y += py;
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
void
centroidC2( const FT &px, const FT &py,
            const FT &qx, const FT &qy,
            const FT &rx, const FT &ry,
            FT &x, FT &y)
{
   x = (px + qx + rx)/FT(3);
   y = (py + qy + ry)/FT(3);
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
   x = (px + qx + rx + sx)/FT(4);
   y = (py + qy + ry + sy)/FT(4);
}

template < class FT >
inline
void
line_from_pointsC2(const FT &px, const FT &py,
                   const FT &qx, const FT &qy,
                   FT &a, FT &b, FT &c) 
{
  a = py - qy;
  b = qx - px;
  // Suggested by Serge Pashkov (psw@rt.kiam.ru) for better numeric stability.
  c = -px*a - py*b;
  // c = px*qy - py*qx;
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
  a = FT(2)*(px - qx);
  b = FT(2)*(py - qy);
  c = CGAL_NTS square(qx) + CGAL_NTS square(qy) -
      CGAL_NTS square(px) - CGAL_NTS square(py);
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
      x = (-b-c)/a + FT(i)*b;
      y = FT(1) - FT(i)*a;
    }
  else
    {
      x = FT(1) + FT(i)*b;
      y = -(a+c)/b - FT(i)*a;
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
  circumcenter_translateC2(qx-px, qy-py, rx-px, ry-py, x, y);
  return CGAL_NTS square(x) + CGAL_NTS square(y);
}

template < class FT >
CGAL_KERNEL_INLINE
FT
squared_distanceC2( const FT &px, const FT &py,
                    const FT &qx, const FT &qy)
{
  return CGAL_NTS square(px-qx) + CGAL_NTS square(py-qy);
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
  return det2x2_by_formula(px-rx,py-ry,qx-rx,qy-ry);
}

CGAL_END_NAMESPACE

#endif // CGAL_CONSTRUCTIONS_KERNEL_FTC2_H
