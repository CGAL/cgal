// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_POLYGON_OFFSET_PRED_FTC2_H
#define CGAL_POLYGON_OFFSET_PRED_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

namespace CGAL {

namespace CGAL_SS_i
{

template <class K>
struct No_cache
{
  typedef boost::optional< Rational< typename K::FT > > Event_time;

  bool is_time_cached(std::size_t)
  {
    return false;
  }

  Event_time time(std::size_t)
  {
    CGAL_error();
    return Event_time();
  }

  void set_time(std::size_t, const Event_time&)
  {}

  void reset(std::size_t)
  {}
};

// Given a triple of oriented straight line segments: (e0,e1,e2) such that
// there exists a distance 'et' for which the offsets lines at 'et' (e0',e1',e2') intersect in a single point;
// returns the relative order of 't' w.r.t 'et'.
// PRECONDITION: There exist a positive distance et for which the offset triple intersect at a single point.
template<class K>
Uncertain<Comparison_result> compare_offset_against_isec_timeC2 ( typename K::FT const& t, intrusive_ptr< Trisegment_2<K> > const& tri )
{
  typedef typename K::FT FT ;

  typedef Rational<FT> Rational ;
  typedef Quotient<FT> Quotient ;

  typedef optional<Rational> Optional_rational ;

  Uncertain<Comparison_result> rResult = Uncertain<Comparison_result>::indeterminate();

  No_cache<K> time_cache;
  Optional_rational et_ = compute_offset_lines_isec_timeC2(tri, time_cache);
  if ( et_ )
  {
    Quotient et = et_->to_quotient();

    CGAL_assertion ( CGAL_NTS certified_is_positive(et) ) ;

    rResult = CGAL_NTS certified_compare( Quotient(t), et);
  }

  return rResult ;

}

} // namespace CGAL_SS_i

} // end namespace CGAL

#endif // CGAL_POLYGON_OFFSET_PRED_FTC2_H //
// EOF //

