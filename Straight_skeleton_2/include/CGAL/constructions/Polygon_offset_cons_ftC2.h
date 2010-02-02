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
// The offsets at the given distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class P_FT, class Segment_2, class Trisegment_2>
boost::optional< boost::tuple< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 
                             , typename Kernel_traits<Trisegment_2>::Kernel::Point_2
                             > 
               >
construct_offset_pointC2 ( P_FT const&                          t
                         , Segment_2 const&                     e0
                         , P_FT const&                          w0
                         , Segment_2 const&                     e1
                         , P_FT const&                          w1
                         , intrusive_ptr< Trisegment_2 > const& tri
                         )
{

  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;
          
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  
  FT x1(0.0), y1(0.0) ;
  FT x2(0.0), y2(0.0) ;

  FT wt0 = t * w0 ;
  FT wt1 = t * w1 ;
  
  CGAL_STSKEL_TRAITS_TRACE("Constructing offset point for t0=" << t0 << " t1=" << t1 << " e0=" << s2str(e0) << " e1=" << s2str(e1) << " tri=" << tri ) ;

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
        FT numX = (wt0 - l0->c()) * l1->b() - (wt1 - l1->c()) * l0->b();
        FT numY = (wt0 - l0->c()) * l1->a() - (wt1 - l1->c()) * l0->a();
    
        x1 = x2 = -numX / den ;
        y1 = y2 =  numY / den ;
        
        ok = CGAL_NTS is_finite(x1) && CGAL_NTS is_finite(y1) ;
      }
      else
      {
        CGAL_STSKEL_TRAITS_TRACE("  DEGENERATE case: Collinear segments involved. Seed event " << ( !tri ? " ABSENT" : " exists." ) ) ;

        Optional_point_2 q = tri ? construct_offset_lines_isecC2(tri) : compute_oriented_midpoint(e0,e1) ;
        
        if ( q )
        {
          CGAL_STSKEL_TRAITS_TRACE("  Seed point: " << p2str(*q) ) ;
           
          FT px, py ;
          line_project_pointC2(l0->a(),l0->b(),l0->c(),q->x(),q->y(),px,py); 
          
          CGAL_STSKEL_TRAITS_TRACE("  Projected seed point: (" << px << "," << py << ")" ) ;

          x1 = px + l0->a() * wt0  ;
          y1 = py + l0->b() * wt0  ;
          
          if ( wt0 != wt1 )
          {
            x2 = px + l0->a() * wt1  ;
            y2 = py + l0->b() * wt1  ;
          }
          else
          {
            x2 = x1 ;
            y2 = y1 ;
          }
          
          ok = CGAL_NTS is_finite(x1) && CGAL_NTS is_finite(y1) && CGAL_NTS is_finite(x2) && CGAL_NTS is_finite(y2);
        }    
      }
    }
  }

  CGAL_STSKEL_TRAITS_TRACE("  RESULT: (" << x << "," << y << ")" << ( ok ? "" : " NONE really" ) ) ;

  return cgal_make_optional(ok, boost::make_tuple( Point_2(x1,y1), Point_2(x2,y2) ) ) ;
}

// Given an offset distance 't' and one oriented line segment e, returns the coordinates (x,y)
// of the source or target point of 'e' offsetted at the given distance.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class P_FT, class Segment_2, class Trisegment_2>
optional< typename Kernel_traits<Trisegment_2>::Kernel::Point_2 > construct_offset_pointC2 ( P_FT const&                          t
                                                                                           , bool                                 at_source 
                                                                                           , Segment_2 const&                     e
                                                                                           , P_FT const&                          w
                                                                                           , intrusive_ptr< Trisegment_2 > const& tri
                                                                                           )
{
  typedef typename Kernel_traits<Trisegment_2>::Kernel K ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;
  typedef typename K::Line_2  Line_2 ;
          
  typedef optional<Point_2> Optional_point_2 ;
  typedef optional<Line_2>  Optional_line_2 ;
  
  FT x(0.0),y(0.0) ;
  
  FT wt = t * w ;

  CGAL_STSKEL_TRAITS_TRACE("Constructing offset point for t=" << t << " e=" << s2str(e) << ( at_source ? " at souce" : " at target" ) << " tri=" << tri ) ;

  Optional_line_2 l = compute_normalized_line_ceoffC2(e) ;

  bool ok = false ;
  
  if ( l )
  {
    CGAL_STSKEL_TRAITS_TRACE("  DEGENERATE case: Collinear segments involved. Seed event " << ( !tri ? " ABSENT" : " exists." ) ) ;

    Optional_point_2 q = tri ? construct_offset_lines_isecC2(tri) : at_source ? e.source() : e.target() ;
    
    if ( q )
    {
      CGAL_STSKEL_TRAITS_TRACE("  Seed point: " << p2str(*q) ) ;
       
      FT px, py ;
      line_project_pointC2(l->a(),l->b(),l->c(),q->x(),q->y(),px,py); 
      
      CGAL_STSKEL_TRAITS_TRACE("  Projected seed point: (" << px << "," << py << ")" ) ;

      x = px + l->a() * wt  ;
      y = py + l->b() * wt  ;
      
      ok = CGAL_NTS is_finite(x) && CGAL_NTS is_finite(y) ;
    }    
  }

  CGAL_STSKEL_TRAITS_TRACE("  RESULT: (" << x << "," << y << ")" << ( ok ? "" : " NONE really" ) ) ;

  return cgal_make_optional(ok, Point_2(x,y)) ;
}

} // namespace CGAL_SS_i

CGAL_END_NAMESPACE

#endif // CGAL_POLYGON_OFFSET_CONS_FTC2_H //
// EOF //
