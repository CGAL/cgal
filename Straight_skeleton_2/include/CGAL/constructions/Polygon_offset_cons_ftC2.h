// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_POLYGON_OFFSET_CONS_FTC2_H
#define CGAL_POLYGON_OFFSET_CONS_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

namespace CGAL {

namespace CGAL_SS_i
{

//TODO: shall we replace it with compute_line_ceoffC2
// Given an oriented 2D straight line segment 'e', computes the normalized coefficients (a,b,c) of the
// supporting line.
// POSTCONDITION: [a,b] is the leftward normal _unit_ (a�+b�=1) vector.
// POSTCONDITION: In case of overflow, an empty optional<> is returned.
template<class K>
optional< Line_2<K> > compute_normalized_line_ceoffC2( Segment_2<K> const& e )
{
  bool finite = true ;
  
  typedef typename K::FT FT ;
  
  FT a (0.0),b (0.0) ,c(0.0)  ;

  if(e.source().y() == e.target().y())
  {
    a = 0 ;
    if(e.target().x() > e.source().x())
    {
      b = 1;
      c = -e.source().y();
    }
    else if(e.target().x() == e.source().x())
    {
      b = 0;
      c = 0;
    }
    else
    {
      b = -1;
      c = e.source().y();
    }

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for HORIZONTAL line:\n" 
                            << s2str(e) 
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else if(e.target().x() == e.source().x())
  {
    b = 0;
    if(e.target().y() > e.source().y())
    {
      a = -1;
      c = e.source().x();
    }
    else if (e.target().y() == e.source().y())
    {
      a = 0;
      c = 0;
    }
    else
    {
      a = 1;
      c = -e.source().x();
    }

    CGAL_STSKEL_TRAITS_TRACE("Line coefficients for VERTICAL line:\n"
                            << s2str(e) 
                            << "\na="<< n2str(a) << ", b=" << n2str(b) << ", c=" << n2str(c)
                           ) ;
  }
  else
  {
    FT sa = e.source().y() - e.target().y();
    FT sb = e.target().x() - e.source().x();
    FT l2 = (sa*sa) + (sb*sb) ;

    if ( CGAL_NTS is_finite(l2) )
    {
      FT l = CGAL_SS_i :: inexact_sqrt(l2);
  
      a = sa / l ;
      b = sb / l ;
  
      c = -e.source().x()*a - e.source().y()*b;
      
      CGAL_STSKEL_TRAITS_TRACE("Line coefficients for line:\n"
                               << s2str(e) 
                               << "\nsa="<< n2str(sa) << "\nsb=" << n2str(sb) << "\nl2=" << n2str(l2) << "\nl=" << n2str(l)
                               << "\na="<< n2str(a) << "\nb=" << n2str(b) << "\nc=" << n2str(c)
                               ) ;
    }
    else finite = false ;
    
  }
  
  if ( finite )
    if ( !CGAL_NTS is_finite(a) || !CGAL_NTS is_finite(b) || !CGAL_NTS is_finite(c) ) 
      finite = false ;

  return cgal_make_optional( finite, K().construct_line_2_object()(a,b,c) ) ;
}

// Given an offset distance 't' and 2 oriented line segments e0 and e1, returns the coordinates (x,y)
// of the intersection of their offsets at the given distance.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a�+b�==1 and (a,b) being the leftward normal vector
// The offsets at the given distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K>
optional< Point_2<K> > construct_offset_pointC2 ( typename K::FT const&                   t
                                                , Segment_2<K> const&                     e0
                                                , Segment_2<K> const&                     e1
                                                , intrusive_ptr< Trisegment_2<K> > const& tri
                                                )
{
  typedef typename K::FT FT ;
  
  typedef Point_2<K>  Point_2 ;
  typedef Line_2<K>   Line_2 ;
          
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  FT x(0.0),y(0.0) ;

  CGAL_STSKEL_TRAITS_TRACE("Constructing offset point for t=" << t << " e0=" << s2str(e0) << " e1=" << s2str(e1) << " tri=" << tri ) ;

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
        CGAL_STSKEL_TRAITS_TRACE("  DEGENERATE case: Collinear segments involved. Seed event " << ( !tri ? " ABSENT" : " exists." ) ) ;

        Optional_point_2 q = tri ? construct_offset_lines_isecC2(tri)->pt : compute_oriented_midpoint(e0,e1) ;
        
        if ( q )
        {
          CGAL_STSKEL_TRAITS_TRACE("  Seed point: " << p2str(*q) ) ;
           
          FT px, py ;
          line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(),px,py); 
          
          CGAL_STSKEL_TRAITS_TRACE("  Projected seed point: (" << px << "," << py << ")" ) ;

          x = px + l0->a() * t  ;
          y = py + l0->b() * t  ;
          
          ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
        }    
      }
    }
  }

  CGAL_STSKEL_TRAITS_TRACE("  RESULT: (" << x << "," << y << ")" << ( ok ? "" : " NONE really" ) ) ;

  return cgal_make_optional(ok,K().construct_point_2_object()(x,y)) ;
}

} // namespace CGAL_SS_i

} // end namespace CGAL

#endif // CGAL_POLYGON_OFFSET_CONS_FTC2_H //
// EOF //

