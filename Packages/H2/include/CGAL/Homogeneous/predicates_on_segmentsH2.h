// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Homogeneous/predicates_on_segmentsH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_PREDICATES_ON_SEGMENTSH2_H
#define CGAL_PREDICATES_ON_SEGMENTSH2_H

#include <CGAL/Homogeneous/SegmentH2.h>
#include <CGAL/Homogeneous/predicates_on_pointsH2.h>

CGAL_BEGIN_NAMESPACE 

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_slopes(const SegmentH2<R>& s1, const SegmentH2<R>& s2)
{
   typedef typename R::FT        FT;

   Comparison_result cmp_y1 = compare_y(s1.source(), s1.target());
   if (cmp_y1 == EQUAL) // horizontal
   {
      Comparison_result cmp_x2 = compare_x(s2.source(), s2.target());

      if (cmp_x2 == EQUAL) return SMALLER;
      FT s_hw = s2.source().hw();
      FT t_hw = s2.target().hw();
      return Comparison_result (
          - CGAL_NTS sign((s2.source().hy()*t_hw - s2.target().hy()*s_hw) *
                          (s2.source().hx()*t_hw - s2.target().hx()*s_hw)) );
   }

   Comparison_result cmp_y2 = compare_y(s2.source(), s2.target());
   if (cmp_y2 == EQUAL)
   {
      Comparison_result cmp_x1 = compare_x(s1.source(), s1.target());

      if (cmp_x1 == EQUAL) return LARGER;
      FT s_hw = s1.source().hw();
      FT t_hw = s1.target().hw();
      return Comparison_result (
             CGAL_NTS sign((s1.source().hy()*t_hw - s1.target().hy()*s_hw) *
                           (s1.source().hx()*t_hw - s1.target().hx()*s_hw)) );
   }

   Comparison_result cmp_x1 = compare_x(s1.source(), s1.target());
   Comparison_result cmp_x2 = compare_x(s2.source(), s2.target());
   if (cmp_x1 == EQUAL)
      return cmp_x2 == EQUAL ? EQUAL : LARGER;

   if (cmp_x2 == EQUAL) return SMALLER;

   FT s1_s_hw = s1.source().hw();
   FT s1_t_hw = s1.target().hw();
   FT s2_s_hw = s2.source().hw();
   FT s2_t_hw = s2.target().hw();
   FT s1_xdiff = s1.source().hx()*s1_t_hw - s1.target().hx()*s1_s_hw;
   FT s1_ydiff = s1.source().hy()*s1_t_hw - s1.target().hy()*s1_s_hw;
   FT s2_xdiff = s2.source().hx()*s2_t_hw - s2.target().hx()*s2_s_hw;
   FT s2_ydiff = s2.source().hy()*s2_t_hw - s2.target().hy()*s2_s_hw;
   Sign s1_sign = CGAL_NTS sign(s1_ydiff * s1_xdiff);
   Sign s2_sign = CGAL_NTS sign(s2_ydiff * s2_xdiff);
 
   if (s1_sign < s2_sign) return SMALLER;
   if (s1_sign > s2_sign) return LARGER;

   if (s1_sign > 0)
     return Comparison_result(
             CGAL_NTS sign ( CGAL_NTS abs(s1_ydiff * s2_xdiff) -
                             CGAL_NTS abs(s2_ydiff * s1_xdiff)) );

   return Comparison_result(
            CGAL_NTS sign ( CGAL_NTS abs(s2_ydiff * s1_xdiff) -
                            CGAL_NTS abs(s1_ydiff * s2_xdiff)) );
}

template < class R >
Comparison_result
compare_y_at_x(const PointH2<R> &p,
               const SegmentH2<R> &s)
{
    // compares the y-coordinates of p and the vertical projection of p on s.
    // Precondition : p is in the x-range of s.

    if (compare_x(s.source(), s.target()) == SMALLER) {
        CGAL_kernel_precondition(compare_x(s.source(), p) != LARGER
		              && compare_x(p, s.target()) != LARGER);
        return (Comparison_result) orientation(p, s.source(), s.target());
    }
    else if (compare_x(s.source(), s.target()) == LARGER) {
        CGAL_kernel_precondition(compare_x(s.target(), p) != LARGER
		              && compare_x(p, s.source()) != LARGER);
        return (Comparison_result) orientation(p, s.target(), s.source());
    }
    else {
        CGAL_kernel_precondition(compare_x(s.target(), p) == EQUAL);
	if (compare_y(p, s.source()) == SMALLER &&
	    compare_y(p, s.target()) == SMALLER)
	    return SMALLER;
	if (compare_y(p, s.source()) == LARGER &&
	    compare_y(p, s.target()) == LARGER)
	    return LARGER;
	return EQUAL;
    }
}

template < class R >
Comparison_result
compare_y_at_x(const PointH2<R> &p,
               const SegmentH2<R> &s1,
               const SegmentH2<R> &s2)
{
    // compares the y-coordinates of the vertical projections of p on s1 and s2
    // Precondition : p is in the x-range of s1 and s2.
    // - if one or two segments are vertical :
    //   - if the segments intersect, return EQUAL
    //   - if not, return the obvious SMALLER/LARGER.

    typedef typename R::FT FT;
    FT px = p.x();
    FT s1sx = s1.source().x();
    FT s1sy = s1.source().y();
    FT s1tx = s1.target().x();
    FT s1ty = s1.target().y();
    FT s2sx = s2.source().x();
    FT s2sy = s2.source().y();
    FT s2tx = s2.target().x();
    FT s2ty = s2.target().y();

    CGAL_kernel_precondition(px >= min(s1sx, s1tx) && px <= max(s1sx, s1tx));
    CGAL_kernel_precondition(px >= min(s2sx, s2tx) && px <= max(s2sx, s2tx));

    if (s1sx != s1tx && s2sx != s2tx) {
	FT s1stx = s1sx-s1tx;
	FT s2stx = s2sx-s2tx;

	return Comparison_result(
	    CGAL_NTS compare(s1sx, s1tx) *
	    CGAL_NTS compare(s2sx, s2tx) *
	    CGAL_NTS compare(-(s1sx-px)*(s1sy-s1ty)*s2stx,
		             (s2sy-s1sy)*s2stx*s1stx
		             -(s2sx-px)*(s2sy-s2ty)*s1stx ));
    }
    else {
	if (s1sx == s1tx) { // s1 is vertical
	    Comparison_result c1, c2;
	    c1 = compare_y_at_x(s1.source(), s2);
	    c2 = compare_y_at_x(s1.target(), s2);
	    if (c1 == c2)
		return c1;
	    return EQUAL;
	}
	// s2 is vertical
	Comparison_result c3, c4;
	c3 = compare_y_at_x(s2.source(), s1);
	c4 = compare_y_at_x(s2.target(), s1);
	if (c3 == c4)
	    return opposite(c3);
	return EQUAL;
    }
}

CGAL_END_NAMESPACE 

#endif // CGAL_PREDICATES_ON_SEGMENTSH2_H
