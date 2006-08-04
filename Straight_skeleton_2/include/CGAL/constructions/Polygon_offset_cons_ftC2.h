// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL$
// $Id$
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_POLYGON_OFFSET_CONS_FTC2_H
#define CGAL_POLYGON_OFFSET_CONS_FTC2_H 1

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i
{

// Given an offset distance 't' and 2 oriented line segments e0 and e1, returns the coordinates (x,y)
// of the intersection of their offsets at the given distance.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at the given distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > construct_offset_pointC2 ( typename K::FT const&         t
                                                , Segment_2<K> const&           e0
                                                , Segment_2<K> const&           e1
                                                , Seeded_trisegment_2<K> const& st
                                                )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K>  Point_2 ;
  typedef Line_2<K>   Line_2 ;
          
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  FT x(0.0),y(0.0) ;

  Optional_line_2 l0 = compute_normalized_line_ceoffC2(e0) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(e1) ;

  bool ok = false ;
  
  if ( l0 && l1 )
  {
    FT den = l1->a() * l0->b() - l0->a() * l1->b() ;
  
    if ( CGAL_NTS is_finite(den) )
    {
      if ( ! CGAL_NTS is_zero(den) )
      {
        FT numX = t * l1->b() - t * l0->b() + l0->b() * l1->c() - l1->b() * l0->c() ;
        FT numY = t * l1->a() - t * l0->a() + l0->a() * l1->c() - l1->a() * l0->c() ;
          
        x = -numX / den ;
        y =  numY / den ;
        
        ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
      }
      else
      {
        Optional_point_2 q = st.event().is_null() ? cgal_make_optional(e1.source()) : construct_offset_lines_isecC2(st);
        if ( q )
        {
          x = q->x() + l0->a() * t  ;
          y = q->y() + l0->b() * t  ;
          
          ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
        }    
      }
      
    }
  }

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

} // namespace CGAL_SS_i

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_CONS_FTC2_H //
// EOF //

