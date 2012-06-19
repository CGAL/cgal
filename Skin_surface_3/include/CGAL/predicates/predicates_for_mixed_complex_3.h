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
// Author(s)     : Nico Kruithof

#ifndef CGAL_PREDICATES_FOR_MIXED_COMPLEX_3_H
#define CGAL_PREDICATES_FOR_MIXED_COMPLEX_3_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

namespace CGAL {

/* q: circumcenter(p1,p2)
 *   <p1-x,p1-q> - s*<p1-q,p1-q>
 * = <p1-x-s(p1-q),p1-q>
 * = <(1-s)p1+s*q-x,p1-q>
 */
template < class FT>
Sign
side_of_mixed_cellC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &p1w,
		     const FT &p2x, const FT &p2y, const FT &p2z, const FT &p2w,
		     const FT &xx, const FT &xy, const FT &xz, 
		     const FT &s) {
  FT x, y, z, dx, dy, dz;
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  x,y,z);
  dx = p1x - p2x;  dy = p1y - p2y;  dz = p1z - p2z;

  return CGAL_NTS sign(((1-s)*p1x+s*x-xx)*(dx) +
		       ((1-s)*p1y+s*y-xy)*(dy) +
		       ((1-s)*p1z+s*z-xz)*(dz));
}

template < class FT>
Sign
side_of_mixed_cellC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &p1w,
		     const FT &p2x, const FT &p2y, const FT &p2z, const FT &p2w,
		     const FT &p3x, const FT &p3y, const FT &p3z, const FT &p3w,
		     const FT &xx, const FT &xy, const FT &xz, 
		     const FT &s) {
  // n is perpendicular to (p1,p2) in the plane of the triangle
  // q is the orthocenter of (p1,p2,p3)
  // t are temporary
  FT nx, ny, nz, qx, qy, qz, tx,ty,tz;

  // t = (p2-p1) x (p3-p1)
  tx = (p2y-p1y) * (p3z-p1z) - (p2z-p1z) * (p3y-p1y);
  ty = (p2z-p1z) * (p3x-p1x) - (p2x-p1x) * (p3z-p1z);
  tz = (p2x-p1x) * (p3y-p1y) - (p2y-p1y) * (p3x-p1x);

  // n = (p2-p1) x t
  nx = (p2y-p1y) * tz - (p2z-p1z) * ty;
  ny = (p2z-p1z) * tx - (p2x-p1x) * tz;
  nz = (p2x-p1x) * ty - (p2y-p1y) * tx;

  // compute q:
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  p3x, p3y, p3z, p3w,
			  qx,qy,qz);

  return CGAL_NTS sign(((1-s)*p1x+s*qx-xx)*nx +
		       ((1-s)*p1y+s*qy-xy)*ny +
		       ((1-s)*p1z+s*qz-xz)*nz);
}

template < class FT>
Sign
side_of_mixed_cellC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &p1w,
		     const FT &p2x, const FT &p2y, const FT &p2z, const FT &p2w,
		     const FT &p3x, const FT &p3y, const FT &p3z, const FT &p3w,
		     const FT &p4x, const FT &p4y, const FT &p4z, const FT &p4w,
		     const FT &xx, const FT &xy, const FT &xz, 
		     const FT &s) {
  // n is perpendicular to (p1,p2,p3) (outward of the tetrahedron)
  // q is the orthocenter of (p1,p2,p3)

  FT nx, ny, nz, qx, qy, qz;

  // n = (p3-p1) x (p2-p1)
  nx = (p3y-p1y) * (p2z-p1z) - (p3z-p1z) * (p2y-p1y);
  ny = (p3z-p1z) * (p2x-p1x) - (p3x-p1x) * (p2z-p1z);
  nz = (p3x-p1x) * (p2y-p1y) - (p3y-p1y) * (p2x-p1x);

  //CGAL_assertion(nx*(p3x-p1x) + ny*(p3y-p1y) + nz*(p3z-p1z) == 0);
  //CGAL_assertion(nx*(p2x-p1x) + ny*(p2y-p1y) + nz*(p2z-p1z) == 0);
  //CGAL_assertion(nx*(p3x-p2x) + ny*(p3y-p2y) + nz*(p3z-p2z) == 0);

  // compute q:
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  p3x, p3y, p3z, p3w,
			  p4x, p4y, p4z, p4w,
			  qx,qy,qz);

  // First term makes up for the orientation of (p1,p2,p3,p4)
  // Second term is the actual sign of test.
  return 
    CGAL_NTS sign(nx*(p1x-p4x) + ny*(p1y-p4y) + nz*(p1z-p4z)) *
    CGAL_NTS sign(((1-s)*p1x+s*qx-xx)*nx +
		  ((1-s)*p1y+s*qy-xy)*ny +
		  ((1-s)*p1z+s*qz-xz)*nz);
}

} //namespace CGAL

#endif // CGAL_PREDICATES_FOR_MIXED_COMPLEX_3_H
