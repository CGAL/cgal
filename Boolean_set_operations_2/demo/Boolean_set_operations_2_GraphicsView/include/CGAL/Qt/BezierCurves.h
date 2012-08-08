// Copyright (c) 2009  GeometryFactory Sarl (France).
// All rights reserved.
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
// 
//
// Author(s) : Fernando Cacciola <fernando.cacciola@geometryfactory.com>

#ifndef CGAL_QT_BEZIER_CURVES_H
#define CGAL_QT_BEZIER_CURVES_H

#include <CGAL/value_type_traits.h>
#include <CGAL/Qt/PiecewiseSetGraphicsItem.h>
#include <CGAL/Qt/BoundaryPiecesGraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#define USE_CLIPPING

namespace CGAL {

namespace Qt {

struct Bezier_helper
{
  template<class OPoint, class IPoint>  static OPoint cvtp2 ( IPoint const& p ) { return OPoint( to_double(p.x()), to_double(p.y()) ) ; }
      
  template<class Bezier_point>
  static Point_2< Simple_cartesian<double> > refine_bezier_point ( Bezier_point const& aP, double aError = 1e-2 )
  {
    typedef typename Bezier_point::Bounding_traits Bounding_traits ;
    typedef typename Bounding_traits::NT           NT ;

    NT min_x, min_y, max_x, max_y ;

    NT const error(aError);
        
    bool lCanRefine = true ;
        
    do
    {
      aP.get_bbox (min_x, min_y, max_x, max_y);

      lCanRefine = !aP.is_exact() && ( CGAL_NTS abs( max_x - min_x) > error || CGAL_NTS abs( max_y - min_y ) > error ) ;
      if ( lCanRefine )
        lCanRefine = aP.refine();
    }
    while ( lCanRefine ) ;
        
    NT const two(2.0);
        
    double x = to_double( ( min_x + max_x ) / two ) ;
    double y = to_double( ( min_y + max_y ) / two ) ;

    return Point_2< Simple_cartesian<double> >(x,y);
  }  

  template<class Bezier_point, class Bezier_curve>
  static typename Bezier_point::Bez_point_bound::NT get_endpoint_parameter( Bezier_point const& aP, Bezier_curve const& aCurve, unsigned int aXID )
  {
    typedef typename Bezier_point::Originator_iterator Originator_iterator ;
    
    Originator_iterator lOrg = aP.get_originator(aCurve, aXID);

    typedef typename Bezier_point::Bez_point_bound::NT NT ;
    
    refine_bezier_point(aP);
        
    return ( lOrg->point_bound().t_min + lOrg->point_bound().t_max ) / NT(2.0);
  }

  template<class Ctrl_points, class FT>
  static void clip_impl ( Ctrl_points& rQ, FT aS, FT aT)
  {
    FT const zero(0.0);
    FT const one (1.0);        
        
    if ( aS >= zero || aT <= one )
    {
      int nc = rQ.size() ;
      int ncc = nc - 1 ;
      
      FT lAlfa = aS ;
      FT lBeta = (aT-aS) / ( one - aS ) ;

      for ( int i = 1 ; i <= ncc ; ++ i )
      {
        for ( int j = 0 ; j < ncc ; ++ j )
          rQ[j] = rQ[j] + lAlfa * ( rQ[j+1] - rQ[j] ) ;
          
        for ( int j = nc - i ; j <= ncc ; ++ j )
          rQ[j] = rQ[j-1] + lBeta * ( rQ[j] - rQ[j-1] ) ;
      }
    }
  }

  template<class Bezier_X_monotone_curve, class Ctrl_points>
  static void clip ( Bezier_X_monotone_curve const& aXMCurve, Ctrl_points& rQ )
  {

    typedef typename Bezier_X_monotone_curve::Curve_2             Bezier_curve ;
    typedef typename Bezier_X_monotone_curve::Bounding_traits::NT BoundNT ;
    typedef typename Bezier_X_monotone_curve::Rat_kernel::Point_2 Rat_point_2 ;
        
    Bezier_curve const& lBC = aXMCurve.supporting_curve();
  
    BoundNT lS = get_endpoint_parameter(aXMCurve.source(), lBC, aXMCurve.xid() );
    BoundNT lT = get_endpoint_parameter(aXMCurve.target(), lBC, aXMCurve.xid() );
    
    bool lFwd = lS <= lT ;
    
    BoundNT lMin = lFwd ? lS : BoundNT(1.0) - lS ;
    BoundNT lMax = lFwd ? lT : BoundNT(1.0) - lT ;
    
    std::copy( lBC.control_points_begin(), lBC.control_points_end(), std::back_inserter(rQ) ) ;
    if ( !lFwd )
      std::reverse(rQ.begin(), rQ.end());
        
    clip_impl(rQ, lMin, lMax ); 
  }


  template<class Control_point_in_iterator, class NT, class Control_point_out_iterator> 
  static void recursive_subdivision( Control_point_in_iterator aBeginCtrlPts
                                   , Control_point_in_iterator aEndCtrlPts
                                   , NT aMin
                                   , NT aMid
                                   , NT aMax
                                   , Control_point_out_iterator aOut
                                   , NT aFlatness
                                   )
  {
    typedef typename value_type_traits<Control_point_in_iterator>::type  Control_point ;
    typedef typename value_type_traits<Control_point_out_iterator>::type Sample_point ;

    Control_point lMinP = point_on_Bezier_curve_2(aBeginCtrlPts, aEndCtrlPts, aMin);
    Control_point lMidP = point_on_Bezier_curve_2(aBeginCtrlPts, aEndCtrlPts, aMid);
    Control_point lMaxP = point_on_Bezier_curve_2(aBeginCtrlPts, aEndCtrlPts, aMax);
    
    NT lB = squared_distance(lMinP, lMaxP);    
    NT lH = squared_distance(lMidP, typename Kernel_traits<Control_point>::Kernel::Line_2(lMinP, lMaxP));    
    
    if ( lH > lB * aFlatness ) 
    {
      recursive_subdivision(aBeginCtrlPts, aEndCtrlPts, aMin, (aMin+aMid)/NT(2.0), aMid, aOut, aFlatness);
      
      *aOut ++ = cvtp2<Sample_point>(lMidP) ;
      
      recursive_subdivision(aBeginCtrlPts, aEndCtrlPts, aMid, (aMid+aMax)/NT(2.0), aMax, aOut, aFlatness);
    }
    else
    {
      *aOut ++ = cvtp2<Sample_point>(lMaxP) ;
    }
    
  }

  template<class Bezier_curve, class NT, class Output_iterator>
  static void sample_curve( Bezier_curve const& aCurve, NT aSourceT, NT aTargetT, Output_iterator aOut, NT aFlatness )
  {
    bool lFwd = aSourceT < aTargetT ;
    
    NT lMinT = lFwd ? aSourceT : NT(1.0) - aSourceT ;
    NT lMaxT = lFwd ? aTargetT : NT(1.0) - aTargetT ;
    
    recursive_subdivision( aCurve.control_points_begin(), aCurve.control_points_end(), lMinT, ( lMinT + lMaxT ) / NT(2.0) , lMaxT, aOut, aFlatness ) ;

  }

  template<class Bezier_curve, class Output_iterator, class NT>
  static void sample_curve( Bezier_curve const& aCurve, Output_iterator aOut, NT aFlatness  ) { sample_curve(aCurve, NT(0.0), NT(1.0), aOut, aFlatness); }

  template<class Bezier_X_monotone_curve, class Output_iterator, class NT>
  static void sample_X_monotone_curve( Bezier_X_monotone_curve const& aXMCurve, bool aFwd, Output_iterator aOut, NT aFlatness  )  
  {
    typedef typename value_type_traits<Output_iterator>::type Output_point ;

    //
    // For a vertical subcurve, the points corresponding to the parameter extremes might not be at the correct vertical positions,
    // so the end points are added explicitely.
    //

    *aOut ++ = refine_bezier_point(aXMCurve.source());

    NT lFromT = get_endpoint_parameter(aFwd ? aXMCurve.source() : aXMCurve.target(), aXMCurve.supporting_curve(), aXMCurve.xid() ) ;
    NT lToT   = get_endpoint_parameter(aFwd ? aXMCurve.target() : aXMCurve.source(), aXMCurve.supporting_curve(), aXMCurve.xid() ) ;
    
    sample_curve(aXMCurve.supporting_curve(), lFromT, lToT, aOut, aFlatness ); 

    *aOut ++ = refine_bezier_point(aXMCurve.target());
  }
  
} ;

struct Bezier_bbox
{
  template<class Bezier_curve>
  Bbox_2 operator()( Bezier_curve const& aBC ) const 
  {
    return CGAL::bbox_2(aBC.control_points_begin(), aBC.control_points_end());
  }
} ;

struct Bezier_X_monotone_bbox
{
  template<class Bezier_X_monotone_curve>
  Bbox_2 operator()( Bezier_X_monotone_curve const& aBXMC ) const 
  {
    typedef typename Bezier_X_monotone_curve::Rat_kernel::Point_2 Rat_point_2 ;
    std::vector<Rat_point_2> lQ ;
    Bezier_helper::clip(aBXMC,lQ);
    return CGAL::bbox_2(lQ.begin(), lQ.end());
  }
} ;

struct Draw_bezier_curve
{
  template<class Bezier_curve, class Path>
  void operator()( Bezier_curve const& aBC, Path& aPath, int aIdx ) const 
  {
    typedef typename Bezier_curve::Bounding_traits::NT BoundNT ;
    typedef typename Bezier_curve::Rat_kernel::Point_2 Rat_point_2 ;
        
    typedef std::vector<Rat_point_2> Rat_point_vector ;
        
    typedef Simple_cartesian<double> Linear_kernel ;
       
    typedef Qt::Converter<Linear_kernel> Converter ;
        
    typedef Point_2<Linear_kernel> Linear_point ;
    
    typedef std::vector<Linear_point> Linear_point_vector ;
    
    if ( aBC.number_of_control_points() == 2 )
    {
      if ( aIdx == 0 )
           aPath.moveTo( Bezier_helper::cvtp2<QPointF>(aBC.control_point(0)) ) ;
      else aPath.lineTo( Bezier_helper::cvtp2<QPointF>(aBC.control_point(0)) ) ; 
      
      aPath.lineTo( Bezier_helper::cvtp2<QPointF>(aBC.control_point(1)) ) ;
    }
    else if ( aBC.number_of_control_points() == 3 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( Bezier_helper::cvtp2<QPointF>(aBC.control_point(0)) ) ;
      else aPath.lineTo ( Bezier_helper::cvtp2<QPointF>(aBC.control_point(0)) ) ; 
      
      aPath.quadTo( Bezier_helper::cvtp2<QPointF>(aBC.control_point(1)), Bezier_helper::cvtp2<QPointF>(aBC.control_point(2)) ) ; 
    }
    else if ( aBC.number_of_control_points() == 4 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( Bezier_helper::cvtp2<QPointF>(aBC.control_point(0)) ) ;
      else aPath.lineTo ( Bezier_helper::cvtp2<QPointF>(aBC.control_point(0)) ) ; 
      
      aPath.cubicTo( Bezier_helper::cvtp2<QPointF>(aBC.control_point(1)), Bezier_helper::cvtp2<QPointF>(aBC.control_point(2)), Bezier_helper::cvtp2<QPointF>(aBC.control_point(3)) ) ; 
    }
    else
    {
      Linear_point_vector lSample ;
      
      Bezier_helper::sample_curve(aBC,std::back_inserter(lSample), BoundNT(1e-4) );
      
      for( typename Linear_point_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
      {
        QPointF lP = Bezier_helper::cvtp2<QPointF>(*it) ;
        
        if ( aIdx == 0 && it == lSample.begin() ) 
             aPath.moveTo(lP) ;
        else aPath.lineTo(lP) ; 
      }
    }
  }
} ;

#ifdef USE_CLIPPING

struct Draw_bezier_X_monotone_curve
{
  template<class Bezier_X_monotone_curve, class Path>
  void operator()( Bezier_X_monotone_curve const& aBXMC, Path& aPath, int aIdx ) const 
  {
    typedef typename Bezier_X_monotone_curve::Curve_2             Bezier_curve ;
    typedef typename Bezier_X_monotone_curve::Bounding_traits::NT BoundNT ;
    typedef typename Bezier_X_monotone_curve::Rat_kernel::Point_2 Rat_point_2 ;
        
    typedef std::vector<Rat_point_2> Rat_point_vector ;

    typedef Simple_cartesian<double> Linear_kernel ;
       
    typedef Qt::Converter<Linear_kernel> Converter ;
        
    typedef Point_2<Linear_kernel> Linear_point ;

    typedef std::vector<Linear_point> Linear_point_vector ;
    
    Rat_point_vector lQ ;
    
    Bezier_helper::clip(aBXMC,lQ);
    
    if ( lQ.size() == 2 )
    {
      if ( aIdx == 0 )
           aPath.moveTo( Bezier_helper::cvtp2<QPointF>(lQ[0]) ) ;
      else aPath.lineTo( Bezier_helper::cvtp2<QPointF>(lQ[0]) ) ; 
      
      aPath.lineTo( Bezier_helper::cvtp2<QPointF>(lQ[1]) ) ;
    }
    else if ( lQ.size() == 3 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( Bezier_helper::cvtp2<QPointF>(lQ[0]) ) ;
      else aPath.lineTo ( Bezier_helper::cvtp2<QPointF>(lQ[0]) ) ; 
      
      aPath.quadTo( Bezier_helper::cvtp2<QPointF>(lQ[1]), Bezier_helper::cvtp2<QPointF>(lQ[2]) ) ; 
    }
    else if ( lQ.size() == 4 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( Bezier_helper::cvtp2<QPointF>(lQ[0]) ) ;
      else aPath.lineTo ( Bezier_helper::cvtp2<QPointF>(lQ[0]) ) ; 
      
      aPath.cubicTo( Bezier_helper::cvtp2<QPointF>(lQ[1]), Bezier_helper::cvtp2<QPointF>(lQ[2]), Bezier_helper::cvtp2<QPointF>(lQ[3]) ) ; 
    }
    else
    {
      
      Linear_point_vector lSample ;
      
      Bezier_helper::sample_X_monotone_curve(aBXMC,true,std::back_inserter(lSample), BoundNT(1e-4) );
      
      for( typename Linear_point_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
      {
        QPointF lP = Bezier_helper::cvtp2<QPointF>(*it) ;
        
        if ( aIdx == 0 && it == lSample.begin() ) 
             aPath.moveTo(lP) ;
        else aPath.lineTo(lP) ; 
      }
    }
  }
} ;

#else

struct Draw_bezier_X_monotone_curve
{
  template<class Bezier_X_monotone_curve, class Path>
  void operator()( Bezier_X_monotone_curve const& aBXMC, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
       
    typedef Qt::Converter<Linear_kernel> Converter ;
        
    typedef Point_2<Linear_kernel> Linear_point ;

    typedef std::vector<Linear_point> Linear_point_vector ;
    
    Converter convert ;
    
    Linear_point_vector lSample ;
    
    Bezier_helper::sample_X_monotone_curve(aBXMC,true,std::back_inserter(lSample) );
    
    for( typename Linear_point_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
    {
      QPointF lP = Bezier_helper::cvtp2<QPointF>(*it) ;
      
      if ( aIdx == 0 && it == lSample.begin() ) 
           aPath.moveTo(lP) ;
      else aPath.lineTo(lP) ; 
    }
  }
} ;

#endif

template<class Bezier_boundary_pieces>
class Bezier_boundary_pieces_graphics_item : public Boundary_pieces_graphics_item<Bezier_boundary_pieces,Draw_bezier_curve,Bezier_bbox>
{
  typedef Boundary_pieces_graphics_item<Bezier_boundary_pieces,Draw_bezier_curve,Bezier_bbox> Base ;
  
public :

  Bezier_boundary_pieces_graphics_item( Bezier_boundary_pieces* aPieces ) : Base(aPieces) {}
} ;

template<class Bezier_boundary>
class Bezier_boundary_graphics_item : public Piecewise_boundary_graphics_item<Bezier_boundary,Draw_bezier_X_monotone_curve,Bezier_X_monotone_bbox>
{
  typedef Piecewise_boundary_graphics_item<Bezier_boundary,Draw_bezier_X_monotone_curve,Bezier_X_monotone_bbox> Base ;
  
public :

  Bezier_boundary_graphics_item( Bezier_boundary* aBoundary ) : Base(aBoundary) {}
} ;

template<class Bezier_region>
class Bezier_region_graphics_item : public Piecewise_region_graphics_item<Bezier_region,Draw_bezier_X_monotone_curve,Bezier_X_monotone_bbox>
{

  typedef Piecewise_region_graphics_item<Bezier_region,Draw_bezier_X_monotone_curve,Bezier_X_monotone_bbox> Base ;
  
public:

  Bezier_region_graphics_item(Bezier_region* aRegion ) : Base(aRegion) {}  
} ;

template<class Bezier_set>
class Bezier_set_graphics_item : public Piecewise_set_graphics_item<Bezier_set,Draw_bezier_X_monotone_curve,Bezier_X_monotone_bbox>
{

  typedef Piecewise_set_graphics_item<Bezier_set,Draw_bezier_X_monotone_curve,Bezier_X_monotone_bbox> Base ;
  
public:

  Bezier_set_graphics_item(Bezier_set* aSet) : Base(aSet) {}
} ;


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_BEZIER_CURVES_H
