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
#ifndef CGAL_POLYGON_OFFSET_CONS_FTC2_H
#define CGAL_POLYGON_OFFSET_CONS_FTC2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

#include <boost/optional/optional.hpp>
#include <boost/intrusive_ptr.hpp>

namespace CGAL {

namespace CGAL_SS_i {

// Given an offset distance 't' and 2 oriented line segments e0 and e1, returns the coordinates (x,y)
// of the intersection of their offsets at the given distance.
//
// PRECONDITIONS:
// The line coefficients must be normalized: a²+b²==1 and (a,b) being the leftward normal vector
// The offsets at the given distance do intersect in a single point.
//
// POSTCONDITION: In case of overflow an empty optional is returned.
//
template<class K, class CoeffCache>
boost::optional< Point_2<K> > construct_offset_pointC2 ( typename K::FT const& t,
                                                         Segment_2_with_ID<K> const& e0,
                                                         Segment_2_with_ID<K> const& e1,
                                                         boost::intrusive_ptr< Trisegment_2<K, Segment_2_with_ID<K> > > const& tri,
                                                         CoeffCache& aCoeff_cache)
{
  typedef typename K::FT FT ;

  typedef Point_2<K>  Point_2 ;
  typedef Line_2<K>   Line_2 ;

  typedef boost::optional<Point_2> Optional_point_2 ;
  typedef boost::optional<Line_2>  Optional_line_2 ;

  FT x(0.0),y(0.0) ;

  CGAL_STSKEL_TRAITS_TRACE("Constructing offset point for t=" << t << " e0=" << s2str(e0) << " e1=" << s2str(e1) << " tri=" << tri ) ;

  Optional_line_2 l0 = compute_normalized_line_ceoffC2(e0, aCoeff_cache) ;
  Optional_line_2 l1 = compute_normalized_line_ceoffC2(e1, aCoeff_cache) ;

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

        Optional_point_2 q = tri ? construct_offset_lines_isecC2(tri, aCoeff_cache)
                                 : compute_oriented_midpoint(e0,e1) ;

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

