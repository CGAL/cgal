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
            CGAL_NTS sign((s2.source().hy()*t_hw - s2.target().hy()*s_hw) *
                          (s2.source().hx()*t_hw - s2.target().hx()*s_hw)) );
   }

   Comparison_result cmp_y2 = compare_y(s2.source(), s2.target());
   if (cmp_y2 == EQUAL)
   {
      Comparison_result cmp_x1 = compare_x(s1.source(), s1.target());

      if (cmp_x1 == EQUAL) return LARGER;
      FT s_hw = s2.source().hw();
      FT t_hw = s2.target().hw();
      return Comparison_result (
             CGAL_NTS sign((s1.source().hy()*s_hw - s1.target().hy()*t_hw) *
                           (s1.source().hx()*s_hw - s1.target().hx()*t_hw)) );
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


CGAL_END_NAMESPACE 

#endif // CGAL_PREDICATES_ON_SEGMENTSH2_H
