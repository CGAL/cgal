// Copyright (c) 2003,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_CONSTRUCTIONS_FTC2_H
#define CGAL_APOLLONIUS_GRAPH_2_CONSTRUCTIONS_FTC2_H 1

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>
#include <CGAL/determinant.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

template < class FT >
inline
void
invert_C2(const FT &x, const FT &y, const FT &wt,
	  FT &new_x, FT &new_y, FT &new_wt)
{
  FT p = CGAL::square(x) + CGAL::square(y) - CGAL::square(wt);

  CGAL_assertion( CGAL::is_positive(p) );

  new_x = x / p;
  new_y = -y / p;
  new_wt = wt / p;
}


template < class FT >
void
w_plane_tangent_line_2(const FT &x1, const FT &y1, const FT &w1,
		       const FT &x2, const FT &y2, const FT &w2,
		       const FT &x3, const FT &y3, const FT &w3,
		       FT       & a, FT       & b, FT       & c)
{
  // we assume that the weight w1 is the smallest among w1, w2, w3.

  FT u2, v2;
  FT u3, v3;
  FT r2, r3;

  invert_C2(x2-x1, y2-y1, w2-w1, u2, v2, r2);
  invert_C2(x3-x1, y3-y1, w3-w1, u3, v3, r3);

  FT Du = u2 - u3;
  FT Dv = v2 - v3;
  FT Dr = r2 - r3;

  FT Duv = determinant(u2, v2, u3, v3);
  FT Dur = determinant(u2, r2, u3, r3);
  FT Dvr = determinant(v2, r2, v3, r3);

  FT D1 = CGAL::square(Du) + CGAL::square(Dv);
  FT D1inv = FT(1) / D1;
  FT sqrtD = CGAL::sqrt(D1 - CGAL::square(Dr));

  a = (Du * Dr - Dv * sqrtD) * D1inv;
  b = (Dv * Dr + Du * sqrtD) * D1inv;
  c = (Du * Dur + Dv * Dvr - Duv * sqrtD) * D1inv;
}


template < class FT >
inline
void
z_plane_circumcircle_2(const FT &x1, const FT &y1, const FT &w1,
		       const FT &x2, const FT &y2, const FT &w2,
		       const FT &x3, const FT &y3, const FT &w3,
		       FT       &cx, FT       &cy, FT      &cwt)
{
  // we assume that the weight w1 is the smallest among w1, w2, w3.

  FT a, b, c;
  w_plane_tangent_line_2(x1, y1, w1, x2, y2, w2,
			 x3, y3, w3, a, b, c);

  cx = -a / (FT(2) * c) + x1;
  cy =  b / (FT(2) * c) + y1;

  // this the only part that is computed at vain when only the center
  // is needed.
#if 0
  FT cwt2 = (CGAL::square(a) + CGAL::square(b)) / 
    (FT(4) * CGAL::square(c));
#endif
  //  cwt = CGAL::sqrt(cwt2) - FT(w1);
  cwt = FT(1) / (FT(2) * c) - FT(w1);
}



template < class FT >
void
ad_circumcircleC2(const FT &x1, const FT &y1, const FT &w1,
		  const FT &x2, const FT &y2, const FT &w2,
		  const FT &x3, const FT &y3, const FT &w3,
		  FT       &cx, FT       &cy, FT      &cwt)
{
  if (CGAL::compare(w2, w1) != LARGER &&
      CGAL::compare(w2, w3) != LARGER) {
    z_plane_circumcircle_2(x2, y2, w2,
			   x3, y3, w3,
			   x1, y1, w1,
			   cx, cy, cwt);
    return;
  } else if (CGAL::compare(w3, w1) != LARGER &&
	     CGAL::compare(w3, w2) != LARGER) {
    z_plane_circumcircle_2(x3, y3, w3,
			   x1, y1, w1,
			   x2, y2, w2,
			   cx, cy, cwt);
    return;
  }
  z_plane_circumcircle_2(x1, y1, w1,
			 x2, y2, w2,
			 x3, y3, w3,
			 cx, cy, cwt);
}

template < class FT >
inline
void
ad_circumcenterC2(const FT &x1, const FT &y1, const FT &w1,
		  const FT &x2, const FT &y2, const FT &w2,
		  const FT &x3, const FT &y3, const FT &w3,
		  FT       &cx, FT       &cy)
{
  FT cwt;
  ad_circumcircleC2(x1, y1, w1,
		    x2, y2, w2,
		    x3, y3, w3,
		    cx, cy, cwt);
}


template < class FT >
void
ad_left_bitangent_lineC2(const FT &x1, const FT &y1, const FT &w1,
			 const FT &x2, const FT &y2, const FT &w2,
			 FT       & a, FT       & b, FT       & c)
{
  FT dx = x1 - x2;
  FT dy = y1 - y2;
  FT dw = w1 - w2;

  FT dxw = determinant(x1, w1, x2, w2);
  FT dyw = determinant(y1, w1, y2, w2);
  FT dxy = determinant(x1, y1, x2, y2);

  FT D1 = CGAL::square(dx) + CGAL::square(dy);
  FT invD1 = FT(1) / D1;
  FT D = D1 - CGAL::square(dw);
  FT sqrtD = CGAL::sqrt(D);

  a = -(dx * dw - dy * sqrtD) * invD1;
  b = -(dy * dw + dx * sqrtD) * invD1;
  c = -(dx * dxw + dy * dyw - dxy * sqrtD) * invD1;
}




} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_CONSTRUCTIONS_FTC2_H
