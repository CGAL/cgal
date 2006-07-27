// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

CGAL_BEGIN_NAMESPACE

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
  FT x, y, z, t;
  t = 1-s;
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  x,y,z);
  return CGAL_NTS sign((t*p1x+s*x-xx)*(p1x-x) +
		       (t*p1y+s*y-xy)*(p1y-y) +
		       (t*p1z+s*z-xz)*(p1z-z));
}

template < class FT>
Sign
side_of_mixed_cellC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &p1w,
		     const FT &p2x, const FT &p2y, const FT &p2z, const FT &p2w,
		     const FT &p3x, const FT &p3y, const FT &p3z, const FT &p3w,
		     const FT &xx, const FT &xy, const FT &xz, 
		     const FT &s) {
  FT x1, y1, z1, x2, y2, z2, t;
  t = 1-s;
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  x1,y1,z1);
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  p3x, p3y, p3z, p3w,
			  x2,y2,z2);
  return CGAL_NTS sign((t*x1+s*x2-xx)*(x1-x2) +
		       (t*y1+s*y2-xy)*(y1-y2) +
		       (t*z1+s*z2-xz)*(z1-z2));
}

template < class FT>
Sign
side_of_mixed_cellC3(const FT &p1x, const FT &p1y, const FT &p1z, const FT &p1w,
		     const FT &p2x, const FT &p2y, const FT &p2z, const FT &p2w,
		     const FT &p3x, const FT &p3y, const FT &p3z, const FT &p3w,
		     const FT &p4x, const FT &p4y, const FT &p4z, const FT &p4w,
		     const FT &xx, const FT &xy, const FT &xz, 
		     const FT &s) {
  FT x1, y1, z1, x2, y2, z2, t;
  t = 1-s;
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  p3x, p3y, p3z, p3w,
			  x1,y1,z1);
  weighted_circumcenterC3(p1x, p1y, p1z, p1w,
			  p2x, p2y, p2z, p2w,
			  p3x, p3y, p3z, p3w,
			  p4x, p4y, p4z, p4w,
			  x2,y2,z2);
  return CGAL_NTS sign((t*x1+s*x2-xx)*(x1-x2) +
		       (t*y1+s*y2-xy)*(y1-y2) +
		       (t*z1+s*z2-xz)*(z1-z2));
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_FOR_MIXED_COMPLEX_3_H
