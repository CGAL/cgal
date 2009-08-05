// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Polyline_simplification_2/demo/Polyline_simplification_2/include/CGAL/Qt/Polyline_simplification_2_graphics_item.h $
// $Id: Polyline_simplification_2_graphics_item.h 48710 2009-04-07 21:41:12Z fcacciola $
// 
//
// Author(s) : Fernando Cacciola <Fernando.Cacciola @geometryfactory.com>

#ifndef CGAL_BEZIER_POLYGONAL_SAMPLER_2_H
#define CGAL_BEZIER_POLYGONAL_SAMPLER_2_H

#include <iterator>
#include <deque>

#include <CGAL/Coercion_traits.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

namespace CGALi {

template<class TargetNT, class SourceNT>
inline TargetNT numeric_cast ( SourceNT const& nt )
{
  typename Coercion_traits<TargetNT,SourceNT>::Cast cast ;
//  return cast(nt);
  return TargetNT( to_double(nt) ) ;
}

template<class NT, class BezierPoint, class BezierCurve>
NT get_approximated_bezier_endpoint_parameter_2( BezierPoint const& p, BezierCurve const& curve, unsigned int xid  )
{
  typedef typename BezierPoint::Originator_iterator Originator_iterator ;
  
  Originator_iterator org = p.get_originator (curve, xid);

  NT threshold(0.0001);
   
  NT t_min = numeric_cast<NT>(org->point_bound().t_min) ;
  NT t_max = numeric_cast<NT>(org->point_bound().t_max) ;
  
  bool can_refine = !p.is_exact();
  
  do
  {
    if ( std::abs(t_max - t_min) <= threshold )
      break ;
    
    if ( can_refine )
    {
      can_refine = p.refine();
      
      t_min = numeric_cast<NT>(org->point_bound().t_min) ;
      t_max = numeric_cast<NT>(org->point_bound().t_max) ;
    }  
  }
  while ( can_refine ) ;
  
  return ( t_min + t_max) / NT(2.0) ;
}



template<class ControlPointContainer, class NT, class OutputPointIterator> 
void bezier_recursive_subdivision_2( ControlPointContainer const& aCtrlPts, NT aMin, NT aMid, NT aMax, OutputPointIterator aOut )
{
  typedef typename ControlPointContainer::value_type            Control_point ;
  typedef typename Kernel_traits<Control_point>::Kernel         K ;
  typedef typename K::Line_2                                    Line ;
  typedef typename value_type_traits<OutputPointIterator>::type Output_point_type ;
  typedef typename Kernel_traits<Output_point_type>::Kernel::FT Output_FT;
  
  Control_point lMinP = point_on_Bezier_curve_2(aCtrlPts.begin(), aCtrlPts.end(), aMin);
  Control_point lMidP = point_on_Bezier_curve_2(aCtrlPts.begin(), aCtrlPts.end(), aMid);
  Control_point lMaxP = point_on_Bezier_curve_2(aCtrlPts.begin(), aCtrlPts.end(), aMax);
  
  NT lB = squared_distance(lMinP, lMaxP);    
  NT lH = squared_distance(lMidP, Line(lMinP, lMaxP));    
  
  if ( lH > lB * NT(0.00001) ) 
  {
    bezier_recursive_subdivision_2(aCtrlPts,aMin,(aMin+aMid)/NT(2.0),aMid,aOut);
    
    *aOut ++ = Output_point_type( numeric_cast<Output_FT>(lMidP.x()), numeric_cast<Output_FT>(lMidP.y()) ) ;
//TRACE("      P:" << Output_point_type( numeric_cast<Output_FT>(lMidP.x()), numeric_cast<Output_FT>(lMidP.y()) ) ) ;
    
    bezier_recursive_subdivision_2(aCtrlPts,aMid,(aMid+aMax)/NT(2.0),aMax,aOut);
  }
  else
  {
    *aOut ++ = Output_point_type( numeric_cast<Output_FT>(lMaxP.x()), numeric_cast<Output_FT>(lMaxP.y()) ) ;
//TRACE("      Q:" << Output_point_type( numeric_cast<Output_FT>(lMaxP.x()), numeric_cast<Output_FT>(lMaxP.y()) ) ) ;
  }
  
}


}

template<class BezierCurve, class NT, class OutputPointIterator>
void sample_bezier_curve_2( BezierCurve const& aBC, NT aSourceT, NT aTargetT, OutputPointIterator aOut )
{
  typedef typename value_type_traits<OutputPointIterator>::type Output_point_type ;
  
  typedef std::deque<Output_point_type> Control_points ;
  
  bool lFwd = aSourceT < aTargetT ;
  
  Control_points lCtrlPoints ;
  int nc = aBC.number_of_control_points() ;
  for ( int i = 0 ; i < nc ; ++ i )
  {
    int j = ( lFwd ? i : nc - i - 1 ) ;
    
    Output_point_type lP ( CGALi::numeric_cast<NT>( aBC.control_point(j).x() )
                         , CGALi::numeric_cast<NT>( aBC.control_point(j).y() )
                         ) ;
    
    lCtrlPoints.push_back (lP);
  }
  
  NT lMinT = lFwd ? aSourceT : NT(1.0) - aSourceT ;
  NT lMaxT = lFwd ? aTargetT : NT(1.0) - aTargetT ;
  
TRACE("    Sampling from " << lMinT << " to " << lMaxT ) ;
  CGALi::bezier_recursive_subdivision_2(lCtrlPoints, lMinT, ( lMinT + lMaxT ) / NT(2.0) , lMaxT, aOut ) ;
}

template<class BezierXMonotoneCurve, class OutputPointIterator>
void sample_bezier_X_monotone_curve_2( BezierXMonotoneCurve const& aBXMC, bool aFwd, OutputPointIterator aOut )
{
  typedef typename BezierXMonotoneCurve::Curve_2                Bezier_curve ;
  typedef typename value_type_traits<OutputPointIterator>::type Output_point_type ;
  typedef typename Kernel_traits<Output_point_type>::Kernel     Output_point_kernel ;
  typedef typename Output_point_kernel::FT                      Output_point_FT ;
  
  Bezier_curve const& lBC = aBXMC.supporting_curve();
  
  Output_point_FT lFromT = CGALi::get_approximated_bezier_endpoint_parameter_2<Output_point_FT>(aFwd ? aBXMC.source() : aBXMC.target(), lBC, aBXMC.xid() ) ;
  Output_point_FT lToT   = CGALi::get_approximated_bezier_endpoint_parameter_2<Output_point_FT>(aFwd ? aBXMC.target() : aBXMC.source(), lBC, aBXMC.xid() ) ;
  
TRACE("    Sampling bezier X monotone curve " << ( aBXMC.is_vertical() ? "[VERTICAL]":"") << ( aBXMC.is_directed_right() ? "[RIGHT]":"[LEFT]") << ":\n    source T=" << lFromT << " P=" << aBXMC.source() << "\n    target T=" << lToT << " P=" << aBXMC.target() );
  sample_bezier_curve_2(lBC, lFromT, lToT, aOut ); 
}

template< class BezierPolygon, class Polygon>
void sample_bezier_polygon_2( BezierPolygon const& aBP, Polygon& rPoly )
{
TRACE("  Sampling bezier polygon. Orientation=" << rPoly.orientation() );
  for( typename BezierPolygon::Curve_const_iterator cit = aBP.curves_begin(); cit != aBP.curves_end(); ++ cit )
    sample_bezier_X_monotone_curve_2(*cit, rPoly.orientation() == COUNTERCLOCKWISE, std::back_inserter(rPoly) );    
}

template<class PolygonWihHoles, class BezierPolygonWithHoles>
void sample_bezier_polygon_with_holes_2( BezierPolygonWithHoles const& aBPWH, PolygonWihHoles& rPolyWithHoles )
{
TRACE("Sampling bezier polygon with holes");
  sample_bezier_polygon_2(aBPWH.outer_boundary(), rPolyWithHoles.outer_boundary() ) ;
  
  for( typename BezierPolygonWithHoles::Hole_const_iterator hit = aBPWH.holes_begin(); hit != aBPWH.holes_end(); ++ hit )
  {
    typename PolygonWihHoles::Polygon_2 lHole ;
    sample_bezier_polygon_2(*hit, lHole );    
    rPolyWithHoles.add_hole(lHole);
  }
TRACE("\n");
}

template<class PolygonWithHoles>
struct Bezier_polygon_with_holes_sampler
{
  typedef PolygonWithHoles result_type ;
  
  template<class BezierPolygonWithHoles>
  result_type operator()( BezierPolygonWithHoles const& aBPWH )
  {
    result_type rR ;
    sample_bezier_polygon_with_holes_2(aBPWH, rR);
    return rR ;  
  }
  
} ;

} // namespace CGAL

#endif // CGAL_BEZIER_POLYGONAL_SAMPLER_2_H
