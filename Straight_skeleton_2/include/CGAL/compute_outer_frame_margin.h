// Copyright (c) 2006-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_COMPUTE_OUTER_FRAME_MARGIN_H
#define CGAL_COMPUTE_OUTER_FRAME_MARGIN_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/Polygon_offset_builder_traits_2.h>
#include <CGAL/Kernel_traits.h>

#include <optional>

#include <algorithm>
#include <iterator>

namespace CGAL {

template<class ForwardPointIterator, class WeightIterator, class Traits>
std::optional< typename Traits::FT > compute_outer_frame_margin ( ForwardPointIterator aBegin
                                                                  , ForwardPointIterator aEnd
                                                                  , WeightIterator       aWBegin
                                                                  , WeightIterator       CGAL_assertion_code(aWEnd)
                                                                  , typename Traits::FT  aOffset
                                                                  , Traits const&        aTraits
                                                                  )
{
  typedef typename Traits::Kernel           Kernel ;
  typedef typename Traits::FT               FT ;
  typedef typename Traits::Point_2          Point_2 ;
  typedef typename Traits::Segment_2        Segment_2 ;
  typedef typename Traits::Trisegment_2_ptr Trisegment_2_ptr ;

  Kernel kernel ;

  typename Kernel::Equal_2                    equal             = kernel.equal_2_object();
  typename Kernel::Collinear_2                collinear         = kernel.collinear_2_object();
  typename Kernel::Compute_squared_distance_2 squared_distance  = kernel.compute_squared_distance_2_object();
  typename Kernel::Construct_segment_2        construct_segment = kernel.construct_segment_2_object();

  typedef std::optional<Point_2> OptionalPoint_2 ;

  CGAL_STSKEL_BUILDER_TRACE(2, "Computing outer frame margin..." );

  FT lMaxSDist(0) ;

  WeightIterator lWIt = aWBegin ;
  ForwardPointIterator lLast = std::prev(aEnd) ;

  bool lOverflow = false ;

  for ( ForwardPointIterator lCurr = aBegin ; lCurr != aEnd ; ++ lCurr, ++ lWIt )
  {
    CGAL_assertion(lWIt != aWEnd);

    ForwardPointIterator lPrev = ( lCurr == aBegin ? lLast  : std::prev  (lCurr) ) ;
    ForwardPointIterator lNext = ( lCurr == lLast  ? aBegin : std::next  (lCurr) ) ;

    if ( !equal(*lPrev,*lCurr) && !equal(*lCurr,*lNext) && !collinear(*lPrev,*lCurr,*lNext) )
    {
      Segment_2 lLEdge = construct_segment(*lPrev,*lCurr);
      Segment_2 lREdge = construct_segment(*lCurr,*lNext);

      WeightIterator lNextWeight = ( lCurr == lLast  ? aWBegin : std::next(lWIt) ) ;

      OptionalPoint_2 lP = aTraits.construct_offset_point_2_object()(aOffset,lLEdge,*lWIt,lREdge,*lNextWeight, Trisegment_2_ptr() );

      if ( !lP )
      {
        lOverflow = true ;
        break ;
      }

      FT lSDist = squared_distance(*lCurr,*lP);

      if (    ! CGAL_NTS is_valid ( lSDist )
           || ! CGAL_NTS is_finite( lSDist )
         )
      {
        lOverflow = true ;
        break ;
      }

      if ( lSDist > lMaxSDist )
        lMaxSDist = lSDist ;
    }
  }

  if ( ! lOverflow )
  {
    FT lDist = CGAL_SS_i::inexact_sqrt(lMaxSDist) ;
    double approx = ceil( to_interval(lDist + ( aOffset * FT(1.05) ) ).second );

    // Add a %5 gap, and ceil to get simpler values
    CGAL_STSKEL_BUILDER_TRACE(4, "outer frame margin: " << approx );
    return std::optional<FT> ( approx ) ;
  }

  return std::nullopt;
}


// `Traits` first is to help overload resolution in the 3-argument version (see below)
template<class Traits, class ForwardPointIterator>
std::optional< typename Traits::FT > compute_outer_frame_margin ( ForwardPointIterator aBegin
                                                                  , ForwardPointIterator aEnd
                                                                  , typename Traits::FT  aOffset
                                                                  , Traits const&        aTraits
                                                                  )
{
  typedef typename Traits::FT FT ;
  std::vector<FT> aUWeights(std::distance(aBegin,aEnd), FT(1)) ;
  return compute_outer_frame_margin(aBegin,aEnd,aUWeights.begin(),aUWeights.end(),aOffset,aTraits) ;
}

template<class ForwardPointIterator, class WeightIterator, class FT>
std::optional<FT> compute_outer_frame_margin(ForwardPointIterator aBegin,
                                               ForwardPointIterator aEnd,
                                               WeightIterator aWBegin,
                                               WeightIterator aWEnd,
                                               const FT aOffset)
{
  typedef typename std::iterator_traits<ForwardPointIterator>::value_type Point_2 ;

  typedef typename Kernel_traits<Point_2>::Kernel K;

  Polygon_offset_builder_traits_2<K> traits ;

  return compute_outer_frame_margin(aBegin,aEnd,aWBegin,aWEnd,aOffset,traits);
}

template<class FT, class ForwardPointIterator>
std::optional<FT> compute_outer_frame_margin(ForwardPointIterator aBegin,
                                               ForwardPointIterator aEnd,
                                               const FT aOffset)
{
  typedef typename std::iterator_traits<ForwardPointIterator>::value_type Point_2 ;

  typedef typename Kernel_traits<Point_2>::Kernel K;
  typedef Polygon_offset_builder_traits_2<K> Builder_traits;
  Builder_traits traits ;

  return compute_outer_frame_margin<Builder_traits>(aBegin,aEnd,aOffset,traits);
}

} // namespace CGAL

#endif // CGAL_COMPUTE_OUTER_FRAME_MARGIN_H
