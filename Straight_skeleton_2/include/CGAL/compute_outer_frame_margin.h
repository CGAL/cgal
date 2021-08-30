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

#include <boost/optional/optional.hpp>

#include <algorithm>
#include <iterator>

namespace CGAL {

template<class ForwardPointIterator, class Traits>
boost::optional< typename Traits::FT > compute_outer_frame_margin ( ForwardPointIterator aBegin
                                                                  , ForwardPointIterator aEnd
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

  typedef boost::optional<Point_2> OptionalPoint_2 ;

  FT lMaxSDist(0) ;

  ForwardPointIterator lLast = std::prev(aEnd) ;

  bool lOverflow = false ;

  for ( ForwardPointIterator lCurr = aBegin ; lCurr != aEnd ; ++ lCurr )
  {
    ForwardPointIterator lPrev = ( lCurr == aBegin ? lLast  : std::prev  (lCurr) ) ;
    ForwardPointIterator lNext = ( lCurr == lLast  ? aBegin : std::next  (lCurr) ) ;

    if ( !equal(*lPrev,*lCurr) && !equal(*lCurr,*lNext) && !collinear(*lPrev,*lCurr,*lNext) )
    {
      Segment_2 lLEdge = construct_segment(*lPrev,*lCurr);
      Segment_2 lREdge = construct_segment(*lCurr,*lNext);

      OptionalPoint_2 lP = Construct_offset_point_2(aTraits)(aOffset,lLEdge,lREdge, Trisegment_2_ptr() );

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

    return boost::optional<FT>( lDist + ( aOffset * FT(1.05) ) ) ; // Add a %5 gap
  }
  else
    return boost::optional<FT>();

}

template<class FT, class ForwardPointIterator>
boost::optional<FT> compute_outer_frame_margin ( ForwardPointIterator aBegin, ForwardPointIterator aEnd,
                                                 FT aOffset )
{
  typedef typename std::iterator_traits<ForwardPointIterator>::value_type Point_2 ;

  typedef typename Kernel_traits<Point_2>::Kernel K;

  Polygon_offset_builder_traits_2<K> traits ;

  return compute_outer_frame_margin(aBegin,aEnd,aOffset,traits);
}

} // end namespace CGAL

#endif // CGAL_COMPUTE_OUTER_FRAME_MARGIN_H //
// EOF //


