#ifndef CGAL_ADDITIVELY_WEIGHTED_DELAUNAY_GRAPH_FTC2_H
#define CGAL_ADDITIVELY_WEIGHTED_DELAUNAY_GRAPH_FTC2_H

#include "../../CGAL/functions_on_signs.h"
#include <CGAL/determinant.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

#define awdg_sltcC2  sign_of_left_tritangent_circle
#define awdg_sdltcC2 sign_of_distance_to_left_tritangent_circle

template < class RT >
Sign
sign_of_left_tritangent_circle(const RT &x1, const RT &y1, const RT &w1,
			       const RT &x2, const RT &y2, const RT &w2,
			       const RT &x3, const RT &y3, const RT &w3)
{
  // we assume that the first circle is the one of smallest radius

  // the result ON_BOUNDARY means that the three circles have a common
  // tangent line; the result ON_BOUNDED_SIDE means that
  // the common tangent circle contains the three original ones and
  // the result ON_UNBOUNDED_SIDE means that the three original
  // circles are on the unbounded side of the original one.

  RT x2s = x2 - x1;
  RT y2s = y2 - y1;
  RT w2s = w2 - w1;
  RT p2s = (CGAL_NTS square(x2s)) + (CGAL_NTS square(y2s))
    - (CGAL_NTS square(w2s));

  RT x3s = x3 - x1;
  RT y3s = y3 - y1;
  RT w3s = w3 - w1;
  RT p3s = (CGAL_NTS square(x3s)) + (CGAL_NTS square(y3s))
    - (CGAL_NTS square(w3s));

  RT dxp = det2x2_by_formula(x2s, p2s, x3s, p3s);
  RT dyp = det2x2_by_formula(y2s, p2s, y3s, p3s);
  RT drp = det2x2_by_formula(w2s, p2s, w3s, p3s);

  RT dxr = det2x2_by_formula(x2s, w2s, x3s, w3s);
  RT dyr = det2x2_by_formula(y2s, w2s, y3s, w3s);
  RT dxy = det2x2_by_formula(x2s, y2s, x3s, y3s);

  RT A = dxp * dxr + dyp * dyr;
  RT B = dxy;
  RT C = (CGAL_NTS square(dxp)) +
    (CGAL_NTS square(dyp)) - (CGAL_NTS square(drp));

  return opposite(sign_a_plus_b_x_sqrt_c(A, B, C));
}

template < class RT >
Sign
sign_of_distance_to_left_tritangent_circle(const RT &x1, const RT &y1,
					   const RT &w1,
					   const RT &x2, const RT &y2,
					   const RT &w2,
					   const RT &x3, const RT &y3,
					   const RT &w3,
					   const RT &x4, const RT &y4,
					   const RT &w4)
{
  // we assume that the first circle is the one of smallest radius

  // the result is the opposite of the sign of the distance of the
  // fourth one from the common tangent circle of the first three

  RT x2s = x2 - x1;
  RT y2s = y2 - y1;
  RT w2s = w2 - w1;
  RT p2s = (CGAL_NTS square(x2s)) + (CGAL_NTS square(y2s))
    - (CGAL_NTS square(w2s));

  RT x3s = x3 - x1;
  RT y3s = y3 - y1;
  RT w3s = w3 - w1;
  RT p3s = (CGAL_NTS square(x3s)) + (CGAL_NTS square(y3s))
    - (CGAL_NTS square(w3s));

  RT x4s = x4 - x1;
  RT y4s = y4 - y1;
  RT w4s = w4 - w1;
  RT p4s = (CGAL_NTS square(x4s)) + (CGAL_NTS square(y4s))
    - (CGAL_NTS square(w4s));

  RT dxp = det2x2_by_formula(x2s, p2s, x3s, p3s);
  RT dyp = det2x2_by_formula(y2s, p2s, y3s, p3s);
  RT drp = det2x2_by_formula(w2s, p2s, w3s, p3s);

  RT dxrp = det3x3_by_formula(x2s, w2s, p2s,
			      x3s, w3s, p3s,
			      x4s, w4s, p4s);

  RT dyrp = det3x3_by_formula(y2s, w2s, p2s,
			      y3s, w3s, p3s,
			      y4s, w4s, p4s);

  RT dxyp = det3x3_by_formula(x2s, y2s, p2s,
			      x3s, y3s, p3s,
			      x4s, y4s, p4s);

  RT A = dxrp * dxp + dyrp * dyp;
  RT B = dxyp;
  RT C = (CGAL_NTS square(dxp)) + (CGAL_NTS square(dyp))
    - (CGAL_NTS square(drp));

  return opposite(sign_a_plus_b_x_sqrt_c(A, B, C));
}


template < class RT >inline
Bounded_side
awdg_testC2(const RT &x1, const RT &y1, const RT &w1,
	    const RT &x2, const RT &y2, const RT &w2,
	    const RT &x3, const RT &y3, const RT &w3,
	    const RT &x4, const RT &y4, const RT &w4)
{
  // this test returns if the fourth weighted point is on the
  // unbounded side of the circle, ccw oriented and which is on the
  // unbounded side, corresponding to the first three weighted points

  // the result ON_BOUNDARY merans that the fourth circle touches the
  // common tangent circle of the first three;
  // the result ON_BOUNDED_SIDE means that the fourth circle
  // intersects the bounded side of the tritangent circle
  // the result ON_UNBOUNDED_SIDE means that the fourth circle does
  // not intersect the bounded side of the tritangent circle

  if ( CGAL_NTS compare(w2, w1) != LARGER &&
       CGAL_NTS compare(w2, w3) != LARGER &&
       CGAL_NTS compare(w2, w4) != LARGER) {
    // we have to perform the test {2, 1, 4, 3}
    return Bounded_side(awdg_sdltcC2(x2, y2, w2,
				     x1, y1, w1,
				     x4, y4, w4,
				     x3, y3, w3));
  } else if ( CGAL_NTS compare(w3, w1) != LARGER &&
	      CGAL_NTS compare(w3, w2) != LARGER &&
	      CGAL_NTS compare(w3, w4) != LARGER) {
    // we have to perform the test {3, 4, 2, 1}
    // and reply negatively
    return Bounded_side(opposite(awdg_sdltcC2(x3, y3, w3,
					      x4, y4, w4,
					      x2, y2, w2,
					      x1, y1, w1)));
  } else if ( CGAL_NTS compare(w4, w1) != LARGER &&
	      CGAL_NTS compare(w4, w2) != LARGER &&
	      CGAL_NTS compare(w4, w3) != LARGER) {
    // we have to perform the test {4, 3, 1, 2}
    // and reply negatively
    return Bounded_side(opposite(awdg_sdltcC2(x4, y4, w4,
					      x3, y3, w3,
					      x1, y1, w1,
					      x2, y2, w2)));
  }

  return Bounded_side(awdg_sdltcC2(x1, y1, w1,
				   x2, y2, w2,
				   x3, y3, w3,
				   x4, y4, w4));
}


template < class RT >
Bounded_side
awdg_testC2(const RT &x1, const RT &y1, const RT &w1,
	    const RT &x2, const RT &y2, const RT &w2,
	    const RT &x3, const RT &y3, const RT &w3)
{
  // this test returns the left tritangent circle, ccw oriented
  // contains the weighted points on its unbounded side

  // the result ON_BOUNDARY means that the three circles have a common
  // tangent line; the result ON_BOUNDED_SIDE means that
  // the common tangent circle contains the three original ones and
  // the result ON_UNBOUNDED_SIDE means that the three original
  // circles are on the unbounded side of the cotangent circle.

  if ( CGAL_NTS compare(w2, w1) != LARGER &&
       CGAL_NTS compare(w2, w3) != LARGER) {
    // we have to perform the test {2, 3, 1}
    return Bounded_side(awdg_sltcC2(x2, y2, w2,
				    x3, y3, w3,
				    x1, y1, w1));
  } else if ( CGAL_NTS compare(w3, w1) != LARGER &&
	      CGAL_NTS compare(w3, w2) != LARGER ) {
    // we have to perform the test {3, 1, 2}
    // and reply negatively
    return Bounded_side(awdg_sltcC2(x3, y3, w3,
				    x1, y1, w1,
				    x2, y2, w2));
  }

  return Bounded_side(awdg_sltcC2(x1, y1, w1,
				  x2, y2, w2,
				  x3, y3, w3));
}


CGAL_END_NAMESPACE

#endif // CGAL_ADDITIVELY_WEIGHTED_DELAUNAY_GRAPH_FTC2_H
