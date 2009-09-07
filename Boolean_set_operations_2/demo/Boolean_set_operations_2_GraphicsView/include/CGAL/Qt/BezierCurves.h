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

#include <CGAL/Qt/PiecewiseSetGraphicsItem.h>
#include <CGAL/Qt/BoundaryPiecesGraphicsItem.h>

namespace CGAL {

namespace Qt {

struct Bezier_helper
{
  template<class Bezier_point, class Bezier_curve>
  static double get_approximated_endpoint_parameter( Bezier_point const& p, Bezier_curve const& curve, unsigned int xid  )
  {
    typedef typename Bezier_point::Originator_iterator Originator_iterator ;
    
    Originator_iterator org = p.get_originator (curve, xid);

    double threshold(0.001);
     
    double t_min = CGAL::to_double(org->point_bound().t_min) ;
    double t_max = CGAL::to_double(org->point_bound().t_max) ;
    
    bool can_refine = !p.is_exact();
    
    do
    {
      if ( std::abs(t_max - t_min) <= threshold )
        break ;
      
      if ( can_refine )
      {
        can_refine = p.refine();
        
        t_min = CGAL::to_double(org->point_bound().t_min) ;
        t_max = CGAL::to_double(org->point_bound().t_max) ;
      }  
    }
    while ( can_refine ) ;
    
    return ( t_min + t_max) / 2.0 ;
  }

  template<class Output_point, class Bezier_curve, class OutputIterator>
  static void get_control_points ( Bezier_curve const& aCurve, bool aFwd, OutputIterator aOut )
  {
    int nc = aCurve.number_of_control_points() ;
    
    for ( int i = aFwd ? 0 : nc - 1 ; aFwd ? i < nc : i >= 0 ; aFwd ? ++ i : -- i )
    {
      *aOut++ = Output_point( CGAL::to_double( aCurve.control_point(i).x() )
                            , CGAL::to_double( aCurve.control_point(i).y() )
                            ) ;
    }
  }  
  
  template<class Output_point, class Bezier_curve, class OutputIterator>
  static void clip ( Bezier_curve const& aCurve, double aS, double aT, bool aFwd, OutputIterator aOut )
  {
    std::vector<Output_point> lQ ;
    
    get_control_points<Output_point>(aCurve, aFwd, std::back_inserter(lQ) ) ;
    
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
  
  template<class Output_point, class Bezier_X_monotone_curve, class OutputIterator>
  static void clip ( Bezier_X_monotone_curve const& aBXMC, OutputIterator aOut )
  {
    typedef typename Bezier_X_monotone_curve::Curve_2 Bezier_curve ;
    
    Bezier_curve const& lBC = aBXMC.supporting_curve();
  
    double lS = get_approximated_endpoint_parameter(aBXMC.source(), lBC, aBXMC.xid());
    double lT = get_approximated_endpoint_parameter(aBXMC.target(), lBC, aBXMC.xid());
    
    bool lFwd = lS <= lT ;
    
    double lMin = lFwd ? lS : 1.0 - lS ;
    double lMax = lFwd ? lT : 1.0 - lT ;
    
    clip<Output_point>(lBC, lMin, lMax, lFwd, aOut ); 
  }
} ;

struct Bezier_bbox
{
  template<class Bezier_curve>
  Bbox_2 operator()( Bezier_curve const& aBC ) const 
  {
    typedef Point_2< Simple_cartesian<double> > Linear_point ;
    std::vector<Linear_point> lQ ;
    Bezier_helper::get_control_points<Linear_point>(aBC,true,std::back_inserter(lQ));
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
    Bezier_helper::clip<Linear_point>(aBXMC,std::back_inserter(lQ));
    return CGAL::bbox_2(lQ.begin(), lQ.end());
  }
} ;

struct Draw_bezier_curve
{
  template<class Bezier_curve, class Path, class Converter>
  void operator()( Bezier_curve const& aBC, Path& aPath, Converter aConvert, int aIdx ) const 
  {
    typedef Point_2< Simple_cartesian<double> > Linear_point ;

    std::vector<Linear_point> lQ ;
    Bezier_helper::get_control_points<Linear_point>(aBC,true,std::back_inserter(lQ));
    if ( lQ.size() == 2 )
    {
      if ( aIdx == 0 )
           aPath.moveTo( aConvert(lQ[0]) ) ;
      else aPath.lineTo( aConvert(lQ[0]) ) ; 
      
      aPath.lineTo( aConvert(lQ[1]) ) ;
    }
    else if ( lQ.size() == 3 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( aConvert(lQ[0]) ) ;
      else aPath.lineTo ( aConvert(lQ[0]) ) ; 
      
      aPath.quadTo( aConvert(lQ[1]), aConvert(lQ[2]) ) ; 
    }
    else if ( lQ.size() == 4 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( aConvert(lQ[0]) ) ;
      else aPath.lineTo ( aConvert(lQ[0]) ) ; 
      
      aPath.cubicTo( aConvert(lQ[1]), aConvert(lQ[2]), aConvert(lQ[3]) ) ; 
    }
    else
    {
      typedef std::pair<double,double> XY ;
      typedef std::vector<XY> XY_vector ;
      XY_vector lSample ;
      aBC.sample(0.0,1.0,100, std::back_inserter(lSample) );
      
      for( typename XY_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
      {
        QPointF lP = aConvert( Linear_point(it->first,it->second) ) ;
        
        if ( aIdx == 0 && it == lSample.begin() ) 
             aPath.moveTo(lP) ;
        else aPath.lineTo(lP) ; 
      }
    }
  }
} ;

struct Draw_bezier_X_monotone_curve
{
  template<class Bezier_X_monotone_curve, class Path, class Converter>
  void operator()( Bezier_X_monotone_curve const& aBXMC, Path& aPath, Converter aConvert, int aIdx ) const 
  {
    typedef Point_2< Simple_cartesian<double> > Linear_point ;

    typedef typename Bezier_X_monotone_curve::Curve_2 Bezier_curve ;
    
    std::vector<Linear_point> lQ ;
    Bezier_helper::clip<Linear_point>(aBXMC,std::back_inserter(lQ));
    if ( lQ.size() == 2 )
    {
      if ( aIdx == 0 )
           aPath.moveTo( aConvert(lQ[0]) ) ;
      else aPath.lineTo( aConvert(lQ[0]) ) ; 
      
      aPath.lineTo( aConvert(lQ[1]) ) ;
    }
    else if ( lQ.size() == 3 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( aConvert(lQ[0]) ) ;
      else aPath.lineTo ( aConvert(lQ[0]) ) ; 
      
      aPath.quadTo( aConvert(lQ[1]), aConvert(lQ[2]) ) ; 
    }
    else if ( lQ.size() == 4 )
    {
      if ( aIdx == 0 )
           aPath.moveTo ( aConvert(lQ[0]) ) ;
      else aPath.lineTo ( aConvert(lQ[0]) ) ; 
      
      aPath.cubicTo( aConvert(lQ[1]), aConvert(lQ[2]), aConvert(lQ[3]) ) ; 
    }
    else
    {
      Bezier_curve const& lBC = aBXMC.supporting_curve();
      
      typedef std::pair<double,double> XY ;
      typedef std::vector<XY> XY_vector ;
      XY_vector lSample ;
      
      double lS = Bezier_helper::get_approximated_endpoint_parameter(aBXMC.source(), lBC, aBXMC.xid());
      double lT = Bezier_helper::get_approximated_endpoint_parameter(aBXMC.target(), lBC, aBXMC.xid());
    
      lBC.sample(lS,lT,100, std::back_inserter(lSample) );
      
      for( typename XY_vector::const_iterator it = lSample.begin() ;  it != lSample.end() ; ++ it )
      {
        QPointF lP = aConvert( Linear_point(it->first,it->second) ) ;
        
        if ( aIdx == 0 && it == lSample.begin() ) 
             aPath.moveTo(lP) ;
        else aPath.lineTo(lP) ; 
      }
    }
  }
} ;

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
