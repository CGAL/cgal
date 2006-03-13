// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_COMPUTE_OUTER_FRAME_MARGIN_H
#define CGAL_COMPUTE_OUTER_FRAME_MARGIN_H

#include <vector>
#include <algorithm>

#include <boost/shared_ptr.hpp>

#include <CGAL/algorithm.h>
#include <CGAL/Polygon_offset_builder_traits_2.h>

CGAL_BEGIN_NAMESPACE

template<class ForwardPointIterator, class Traits>
typename Traits::FT compute_outer_frame_margin ( ForwardPointIterator aBegin
                                               , ForwardPointIterator aEnd
                                               , typename Traits::FT  aOffset
                                               , Traits const&        aTraits
                                               )
{
  typedef typename Traits::Kernel  Kernel ;
  typedef typename Traits::Point_2 Point_2 ;
  typedef typename Traits::FT      FT ;
  
  typedef typename Traits::Construct_ss_edge_2::Edge Edge ;
  
  Kernel kernel ;
  
  typename Kernel::Equal_2                    equal            = kernel.equal_2_object();
  typename Kernel::Collinear_2                collinear        = kernel.collinear_2_object();
  typename Kernel::Compute_squared_distance_2 squared_distance = kernel.compute_squared_distance_2_object();
  
  FT lMaxSDist(0.0) ;
  
  ForwardPointIterator lLast = predecessor(aEnd) ;
  
  for ( ForwardPointIterator lCurr = aBegin ; lCurr < aEnd ; ++ lCurr )
  {
    ForwardPointIterator lPrev = ( lCurr == aBegin ? lLast  : predecessor(lCurr) ) ;
    ForwardPointIterator lNext = ( lCurr == lLast  ? aBegin : successor  (lCurr) ) ;
    
    if ( !equal(*lPrev,*lCurr) && !equal(*lCurr,*lNext) && !collinear(*lPrev,*lCurr,*lNext) )
    {
      Edge    lLEdge = Construct_ss_edge_2     <Traits>(aTraits)()(*lPrev,*lCurr);
      Edge    lREdge = Construct_ss_edge_2     <Traits>(aTraits)()(*lCurr,*lNext);
      Point_2 lP     = Construct_offset_point_2<Traits>(aTraits)()(aOffset,lLEdge,lREdge);
     
      FT lSDist = squared_distance(*lCurr,lP);
      
      if ( lSDist > lMaxSDist )
        lMaxSDist = lSDist ;
    }
  }
  
  FT lDist = CGAL_NTS sqrt(lMaxSDist) ;
  
  return lDist + ( aOffset * FT(1.05) ) ; // Add a %5 gap
}                              

template<class ForwardPointIterator, class FT>
FT compute_outer_frame_margin ( ForwardPointIterator aBegin, ForwardPointIterator aEnd, FT aOffset )
{
  typedef typename std::iterator_traits<ForwardPointIterator>::value_type Point_2 ;
  
  typedef typename Kernel_traits<Point_2>::Kernel K; 
  
  Polygon_offset_builder_traits_2<K> traits ;
  
  return compute_outer_frame_margin(aBegin,aEnd,aOffset,traits);
}                                                                 

CGAL_END_NAMESPACE

#endif // CGAL_COMPUTE_EXTERIOR_OFFSET_CONTOUR_FRAME_MARGIN_2_H //
// EOF //

 
