// Copyright (c) 2009  GeometryFactory Sarl (France).
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
  template<class Bezier_point, class Bezier_curve>
  static double get_approximated_endpoint_parameter( Bezier_point const& aP, Bezier_curve const& aCurve, unsigned int aXID, double aError = 1e-3 )
  {
    typedef typename Bezier_point::Originator_iterator Originator_iterator ;
    
    Originator_iterator lOrg = aP.get_originator(aCurve, aXID);

    double t_min = CGAL::to_double(lOrg->point_bound().t_min) ;
    double t_max = CGAL::to_double(lOrg->point_bound().t_max) ;
    
    bool lCanRefine = !aP.is_exact();
    
    do
    {
      if ( std::abs(t_max - t_min) <= aError )
        break ;
      
      if ( lCanRefine )
      {
        lCanRefine = aP.refine();
        
        t_min = CGAL::to_double(lOrg->point_bound().t_min) ;
        t_max = CGAL::to_double(lOrg->point_bound().t_max) ;
      }  
    }
    while ( lCanRefine ) ;
    
    return ( t_min + t_max) / 2.0 ;
  }

  template<class Bezier_curve, class Output_iterator>
  static void get_approximated_control_points ( Bezier_curve const& aCurve, bool aFwd, Output_iterator aOut )
  {
    typedef typename value_type_traits<Output_iterator>::type Output_point ;

    int nc = aCurve.number_of_control_points() ;
    
    for ( int i = aFwd ? 0 : nc - 1 ; aFwd ? i < nc : i >= 0 ; aFwd ? ++ i : -- i )
    {
      *aOut++ = Output_point( CGAL::to_double( aCurve.control_point(i).x() )
                            , CGAL::to_double( aCurve.control_point(i).y() )
                            ) ;
    }
  }  

  template<class Bezier_point>
  static Point_2< Simple_cartesian<double> > get_approximated_bezier_point ( Bezier_point const& aP, double aError = 1e-3 )
  {
    bool lCanRefine = !aP.is_exact();

    typedef typename Bezier_point::Bounding_traits Bounding_traits ;
    typedef typename Bounding_traits::NT           NT ;

    NT min_x, min_y, max_x, max_y ;

    do
    {
      aP.get_bbox (min_x, min_y, max_x, max_y);

      if ( std::abs( to_double(max_x - min_x) ) <= aError && std::abs( to_double(max_y - min_y) ) <= aError )
        break ;

      aP.refine();

    }
    while ( lCanRefine ) ;

    double x = to_double( min_x + max_x ) / 2.0 ;
    double y = to_double( min_y + max_y ) / 2.0 ;

    return Point_2< Simple_cartesian<double> >(x,y);
  }  

  template<class Bezier_curve, class Output_iterator>
  static void approximated_clip ( Bezier_curve const& aCurve, double aS, double aT, bool aFwd, Output_iterator aOut )
  {
    typedef typename value_type_traits<Output_iterator>::type Output_point ;

    std::vector<Output_point> lQ ;
    
    get_approximated_control_points(aCurve, aFwd, std::back_inserter(lQ) ) ;
    
    if ( aS >= 0.0 || aT <= 1.0 )
    {
      int nc = aCurve.number_of_control_points() ;
      int ncc = nc - 1 ;
      
      double lAlfa = aS ;
      double lBeta = (aT-aS) / ( 1.0 - aS ) ;

      for ( int i = 1 ; i <= ncc ; ++ i )
      {
        for ( int j = 0 ; j < ncc ; ++ j )
          lQ[j] = lQ[j] + lAlfa * ( lQ[j+1] - lQ[j] ) ;
          
        for ( int j = nc - i ; j <= ncc ; ++ j )
          lQ[j] = lQ[j-1] + lBeta * ( lQ[j] - lQ[j-1] ) ;
      }
    }
    
    std::copy(lQ.begin(), lQ.end(), aOut );    
    
  }


  template<class Bezier_X_monotone_curve, class Output_iterator>
  static void approximated_clip ( Bezier_X_monotone_curve const& aXMCurve, Output_iterator aOut, double aApproxError = 1e-3 )
  {
    typedef typename value_type_traits<Output_iterator>::type Output_point ;

    typedef typename Bezier_X_monotone_curve::Curve_2 Bezier_curve ;
    
    Bezier_curve const& lBC = aXMCurve.supporting_curve();
  
    double lS = get_approximated_endpoint_parameter(aXMCurve.source(), lBC, aXMCurve.xid(), aApproxError);
    double lT = get_approximated_endpoint_parameter(aXMCurve.target(), lBC, aXMCurve.xid(), aApproxError);
    
    bool lFwd = lS <= lT ;
    
    double lMin = lFwd ? lS : 1.0 - lS ;
    double lMax = lFwd ? lT : 1.0 - lT ;
    
    approximated_clip(lBC, lMin, lMax, lFwd, aOut ); 
  }
  
  template<class Control_point_in_iterator, class NT, class Control_point_out_iterator> 
  static void approximated_recursive_subdivision( Control_point_in_iterator aBeginCtrlPts
                                                , Control_point_in_iterator aEndCtrlPts
                                                , NT aMin
                                                , NT aMid
                                                , NT aMax
                                                , Control_point_out_iterator aOut
                                                , NT aFlatness = 1e-4
                                                )
  {
    typedef typename value_type_traits<Control_point_in_iterator>::type Control_point ;

    Control_point lMinP = point_on_Bezier_curve_2(aBeginCtrlPts, aEndCtrlPts, aMin);
    Control_point lMidP = point_on_Bezier_curve_2(aBeginCtrlPts, aEndCtrlPts, aMid);
    Control_point lMaxP = point_on_Bezier_curve_2(aBeginCtrlPts, aEndCtrlPts, aMax);
    
    NT lB = squared_distance(lMinP, lMaxP);    
    NT lH = squared_distance(lMidP, typename Kernel_traits<Control_point>::Kernel::Line_2(lMinP, lMaxP));    
    
    if ( lH > lB * aFlatness ) 
    {
      approximated_recursive_subdivision(aBeginCtrlPts,aEndCtrlPts,aMin,(aMin+aMid)/NT(2.0),aMid,aOut,aFlatness);
      
      *aOut ++ = lMidP ;
      
      approximated_recursive_subdivision(aBeginCtrlPts,aEndCtrlPts,aMid,(aMid+aMax)/NT(2.0),aMax,aOut,aFlatness);
    }
    else
    {
      *aOut ++ = lMaxP ;
    }
    
  }

  template<class Bezier_curve, class NT, class Output_iterator>
  static void sample_curve( Bezier_curve const& aCurve, NT aSourceT, NT aTargetT, Output_iterator aOut, double aFlatness = 1e-4 )
  {
    typedef typename value_type_traits<Output_iterator>::type Output_point ;
    
    typedef std::vector<Output_point> Control_points ;
    
    Control_points lControlPoints ;
    
    bool lFwd = aSourceT < aTargetT ;
    
    get_approximated_control_points(aCurve, lFwd, std::back_inserter(lControlPoints) ) ;
    
    NT lMinT = lFwd ? aSourceT : NT(1.0) - aSourceT ;
    NT lMaxT = lFwd ? aTargetT : NT(1.0) - aTargetT ;
    
    approximated_recursive_subdivision(lControlPoints.begin(), lControlPoints.end(), lMinT, ( lMinT + lMaxT ) / NT(2.0) , lMaxT, aOut, aFlatness ) ;

  }

  template<class Bezier_curve, class Output_iterator>
  static void sample_curve( Bezier_curve const& aCurve, Output_iterator aOut, double aFlatness = 1e-4 ) { sample_curve(aCurve, 0.0, 1.0, aOut, aFlatness); }

  template<class Bezier_X_monotone_curve, class Output_iterator>
  static void sample_X_monotone_curve( Bezier_X_monotone_curve const& aXMCurve, bool aFwd, Output_iterator aOut, double aFlatness = 1e-4, double aEPApproxError = 1e-3 )  
  {
    typedef typename value_type_traits<Output_iterator>::type Output_point ;

    //
    // For a vertical subcurve, the points corresponding to the parameter extremes might not be at the correct vertical positions,
    // so the end points are added explicitely.
    //

    *aOut ++ = get_approximated_bezier_point(aXMCurve.source());

    double lFromT = get_approximated_endpoint_parameter(aFwd ? aXMCurve.source() : aXMCurve.target(), aXMCurve.supporting_curve(), aXMCurve.xid(), aEPApproxError ) ;
    double lToT   = get_approximated_endpoint_parameter(aFwd ? aXMCurve.target() : aXMCurve.source(), aXMCurve.supporting_curve(), aXMCurve.xid(), aEPApproxError ) ;
    
    sample_curve(aXMCurve.supporting_curve(), lFromT, lToT, aOut, aFlatness ); 

    *aOut ++ = get_approximated_bezier_point(aXMCurve.target());
  }
  
} ;

struct Bezier_bbox
{
  template<class Bezier_curve>
  Bbox_2 operator()( Bezier_curve const& aBC ) const 
  {
    typedef Point_2< Simple_cartesian<double> > Linear_point ;
    std::vector<Linear_point> lQ ;
    Bezier_helper::get_approximated_control_points(aBC,true,std::back_inserter(lQ));
    return CGAL::bbox_2(lQ.begin(), lQ.end());
  }
} ;

struct Bezier_X_monotone_bbox
{
  template<class Bezier_X_monotone_curve>
  Bbox_2 operator()( Bezier_X_monotone_curve const& aBXMC ) const 
  {
    typedef Point_2< Simple_cartesian<double> > Linear_point ;
    std::vector<Linear_point> lQ ;
    Bezier_helper::approximated_clip(aBXMC,std::back_inserter(lQ));
    return CGAL::bbox_2(lQ.begin(), lQ.end());
  }
} ;

struct Draw_bezier_curve
{
  template<class Bezier_curve, class Path>
  void operator()( Bezier_curve const& aBC, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
       
    typedef Qt::Converter<Linear_kernel> Converter ;
        
    typedef Point_2<Linear_kernel> Linear_point ;
    
    typedef std::vector<Linear_point> Linear_point_vector ;
    
    Linear_point_vector lQ ;
    
    Bezier_helper::get_approximated_control_points(aBC,true,std::back_inserter(lQ));
    
    Converter convert ;
    
    if ( lQ.size() == 2 )
    {
      if ( aIdx == 0 )
           aPath.moveTo( convert(lQ[0]) ) ;
      else aPath.lineTo( convert(lQ[0]) ) ; 
      
      aPath.lineTo( convert(lQ[1]) ) ;
    }
    else if ( lQ.size() == 3 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( convert(lQ[0]) ) ;
      else aPath.lineTo ( convert(lQ[0]) ) ; 
      
      aPath.quadTo( convert(lQ[1]), convert(lQ[2]) ) ; 
    }
    else if ( lQ.size() == 4 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( convert(lQ[0]) ) ;
      else aPath.lineTo ( convert(lQ[0]) ) ; 
      
      aPath.cubicTo( convert(lQ[1]), convert(lQ[2]), convert(lQ[3]) ) ; 
    }
    else
    {
      Linear_point_vector lSample ;
      
      Bezier_helper::sample_curve(aBC,std::back_inserter(lSample) );
      
      for( typename Linear_point_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
      {
        QPointF lP = convert(*it) ;
        
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
    typedef Simple_cartesian<double> Linear_kernel ;
       
    typedef Qt::Converter<Linear_kernel> Converter ;
        
    typedef Point_2<Linear_kernel> Linear_point ;

    typedef std::vector<Linear_point> Linear_point_vector ;
    
    typedef typename Bezier_X_monotone_curve::Curve_2 Bezier_curve ;
    
    Linear_point_vector lQ ;

    Converter convert ;
    
    Bezier_helper::approximated_clip(aBXMC,std::back_inserter(lQ));
    
    if ( lQ.size() == 2 )
    {
      if ( aIdx == 0 )
           aPath.moveTo( convert(lQ[0]) ) ;
      else aPath.lineTo( convert(lQ[0]) ) ; 
      
      aPath.lineTo( convert(lQ[1]) ) ;
    }
    else if ( lQ.size() == 3 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( convert(lQ[0]) ) ;
      else aPath.lineTo ( convert(lQ[0]) ) ; 
      
      aPath.quadTo( convert(lQ[1]), convert(lQ[2]) ) ; 
    }
    else if ( lQ.size() == 4 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( convert(lQ[0]) ) ;
      else aPath.lineTo ( convert(lQ[0]) ) ; 
      
      aPath.cubicTo( convert(lQ[1]), convert(lQ[2]), convert(lQ[3]) ) ; 
    }
    else
    {
      
      Linear_point_vector lSample ;
      
      Bezier_helper::sample_X_monotone_curve(aBXMC,true,std::back_inserter(lSample) );
      
      for( typename Linear_point_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
      {
        QPointF lP = convert(*it) ;
        
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
      QPointF lP = convert(*it) ;
      
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
