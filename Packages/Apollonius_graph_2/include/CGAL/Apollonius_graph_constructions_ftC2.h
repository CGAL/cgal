// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/Apollonius_graph_constructions_ftC2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_CONSTRUCTIONS_FTC2_H
#define CGAL_APOLLONIUS_GRAPH_CONSTRUCTIONS_FTC2_H

#include <CGAL/determinant.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
inline
void
invert_C2(const FT &x, const FT &y, const FT &wt,
	  FT &new_x, FT &new_y, FT &new_wt)
{
  FT p = CGAL_NTS square(x) + CGAL_NTS square(y) - CGAL_NTS square(wt);

  CGAL_assertion( CGAL_NTS is_positive(p) );

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

  FT Duv = det2x2_by_formula(u2, v2, u3, v3);
  FT Dur = det2x2_by_formula(u2, r2, u3, r3);
  FT Dvr = det2x2_by_formula(v2, r2, v3, r3);

  FT D1 = CGAL_NTS square(Du) + CGAL_NTS square(Dv);
  FT D1inv = FT(1) / D1;
  FT sqrtD = CGAL_NTS sqrt(D1 - CGAL_NTS square(Dr));

  a = (Du * Dr - Dv * sqrtD) * D1inv;
  b = (Dv * Dr + Du * sqrtD) * D1inv;
  c = (Du * Dur + Dv * Dvr - Duv * sqrtD) * D1inv;
}


template < class FT, class We >
inline
void
z_plane_circumcircle_2(const FT &x1, const FT &y1, const We &w1,
		       const FT &x2, const FT &y2, const We &w2,
		       const FT &x3, const FT &y3, const We &w3,
		       FT       &cx, FT       &cy, We      &cwt)
{
  // we assume that the weight w1 is the smallest among w1, w2, w3.

  FT a, b, c;
  w_plane_tangent_line_2(x1, y1, FT(w1), x2, y2, FT(w2),
			 x3, y3, FT(w3), a, b, c);

  cx = -a / (FT(2) * c) + x1;
  cy =  b / (FT(2) * c) + y1;

  // this the only part that is computed at vain when only the center
  // is needed.
#if 0
  FT cwt2 = (CGAL_NTS square(a) + CGAL_NTS square(b)) / 
    (FT(4) * CGAL_NTS square(c));
#endif
  //  cwt = We(  CGAL_NTS sqrt(cwt2) - FT(w1)  );
  cwt = We( FT(1) / (FT(2) * c) - FT(w1));
}



template < class FT, class We >
inline
void
ad_circumcenterC2(const FT &x1, const FT &y1, const We &w1,
		  const FT &x2, const FT &y2, const We &w2,
		  const FT &x3, const FT &y3, const We &w3,
		  FT       &cx, FT       &cy)
{
  We cwt;
  ad_circumcircleC2(x1, y1, w1,
		    x2, y2, w2,
		    x3, y3, w3,
		    cx, cy, cwt);
}


template < class FT, class We >
void
ad_circumcircleC2(const FT &x1, const FT &y1, const We &w1,
		  const FT &x2, const FT &y2, const We &w2,
		  const FT &x3, const FT &y3, const We &w3,
		  FT       &cx, FT       &cy, We      &cwt)
{
  if (CGAL_NTS compare(w2, w1) != LARGER &&
      CGAL_NTS compare(w2, w3) != LARGER) {
    z_plane_circumcircle_2(x2, y2, w2,
			   x3, y3, w3,
			   x1, y1, w1,
			   cx, cy, cwt);
    return;
  } else if (CGAL_NTS compare(w3, w1) != LARGER &&
	     CGAL_NTS compare(w3, w2) != LARGER) {
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

template < class FT, class We >
void
ad_left_bitangent_lineC2(const FT &x1, const FT &y1, const We &w1,
			 const FT &x2, const FT &y2, const We &w2,
			 FT       & a, FT       & b, FT       & c)
{
  FT dx = x1 - x2;
  FT dy = y1 - y2;
  FT dw = FT(w1 - w2);

  FT dxw = det2x2_by_formula(x1, w1, x2, w2);
  FT dyw = det2x2_by_formula(y1, w1, y2, w2);
  FT dxy = det2x2_by_formula(x1, y1, x2, y2);

  FT D1 = CGAL_NTS square(dx) + CGAL_NTS square(dy);
  FT invD1 = FT(1) / D1;
  FT D = D1 - CGAL_NTS square(dw);
  FT sqrtD = CGAL_NTS sqrt(D);

  a = -(dx * dw - dy * sqrtD) * invD1;
  b = -(dy * dw + dx * sqrtD) * invD1;
  c = -(dx * dxw + dy * dyw - dxy * sqrtD) * invD1;
}






CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_CONSTRUCTIONS_FTC2_H
