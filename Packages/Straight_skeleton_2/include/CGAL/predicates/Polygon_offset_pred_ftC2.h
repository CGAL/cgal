// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/predicates/Polygon_offset_pred_ftC2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_PRED_FTC2_H
#define CGAL_POLYGON_OFFSET_PRED_FTC2_H 1

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SLS_i
{

// Given a triple of oriented lines in _normalized_ implicit form: (e0,e1,e2) such that
// there exists a distance 'et' for which the offsets lines at 'et' (e0',e1',e2') intersect in a single point;
// returns the relative order of t w.r.t et.
// PRECONDITION: There exist a positive distance et for which the offset triple intersect at a single point.
template<class FT>
Uncertain<Comparison_result>
compare_offset_against_isec_timeC2 ( FT t, Triedge<FT> const& triedge )
{
  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  SortedTriedge<FT> sorted = collinear_sort(triedge);

  if ( sorted.is_valid() )
  {
    Quotient<FT> et = compute_offset_lines_isec_timeC2(sorted).to_quotient();

    CGAL_assertion ( CGAL_NTS certified_is_positive(et) ) ;

    rResult = CGAL_NTS certified_compare(et,t);
  }

  return rResult ;

}

} // namespace CGAL_SLS_i

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_PRED_FTC2_H //
// EOF //

