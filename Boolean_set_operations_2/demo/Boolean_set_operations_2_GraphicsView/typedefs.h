// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Boolean_set_operations_2/demo/Boolean_set_operations_2/typedefs.h $
// $Id: typedefs.h 37003 2007-03-10 16:55:12Z spion $
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//
#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

//
// Linear polygons
//
typedef CGAL::Simple_cartesian<double>            Linear_kernel ;
typedef CGAL::Polygon_2<Linear_kernel>            Linear_polygon;
typedef CGAL::Polygon_with_holes_2<Linear_kernel> Linear_polygon_with_holes;

typedef Linear_kernel::Point_2 Linear_point ;

//
// Circlular polygons
//

#ifdef CGAL_USE_GMP

  typedef CGAL::Gmpq                     Base_nt;

#else

  typedef CGAL::Quotient<CGAL::MP_Float> Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt> Coord_type;


struct Circular_kernel : public CGAL::Cartesian<Coord_type> {};

typedef CGAL::Gps_circle_segment_traits_2<Circular_kernel>   Circular_traits;
typedef Circular_traits::Curve_2                             Circular_curve;
typedef Circular_traits::X_monotone_curve_2                  Circular_X_monotone_curve;
typedef Circular_traits::Point_2                             Circular_point ;
typedef Circular_traits::Polygon_2                           Circular_polygon;
typedef CGAL::General_polygon_with_holes_2<Circular_polygon> Circular_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Circular_traits>         Circular_polygon_set;

struct Compute_circular_X_monotone_cuve_bbox
{
  CGAL::Bbox_2 operator()( Circular_X_monotone_curve const& curve ) const 
  {
    return curve.bbox();   
  }
} ;

struct Draw_circular_X_monotone_cuve
{
  template<class Path, class Converter>
  void operator()( Circular_X_monotone_curve const& curve, Path* aPath, Converter aConvert, int aIdx ) const 
  {
    Linear_point lS( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
    Linear_point lT( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
    
    if ( aIdx == 0 )
         aPath->moveTo( aConvert( lS ) ) ;
    else aPath->lineTo( aConvert( lS ) ) ;
    
    aPath->lineTo( aConvert( lT ) ) ;
  }
} ;
typedef CGAL::Qt::GeneralPolygonSetGraphicsItem<Circular_polygon_set,Compute_circular_X_monotone_cuve_bbox,Draw_circular_X_monotone_cuve> Circular_GI;



//
// Bezier curves typedefs
//
#ifdef CGAL_USE_CORE

typedef CGAL::CORE_algebraic_number_traits            Bezier_nt_traits;
typedef Bezier_nt_traits::Rational                    Bezier_rational;
typedef Bezier_nt_traits::Algebraic                   Bezier_algebraic;

struct Bezier_rat_kernel  : public CGAL::Cartesian<Bezier_rational>  {};
struct Bezier_alg_kernel  : public CGAL::Cartesian<Bezier_algebraic> {};

struct Bezier_traits : public CGAL::Arr_Bezier_curve_traits_2<Bezier_rat_kernel, Bezier_alg_kernel, Bezier_nt_traits> {};
  
typedef Bezier_rat_kernel::Point_2                      Bezier_rat_point;
typedef Bezier_traits::Curve_2                          Bezier_curve;
typedef Bezier_traits::X_monotone_curve_2               Bezier_X_monotone_curve;
typedef Bezier_traits::Point_2                          Bezier_point;
typedef CGAL::Gps_traits_2<Bezier_traits>               Bezier_gps_traits;
typedef Bezier_gps_traits::General_polygon_2            Bezier_polygon;
typedef std::vector<Bezier_polygon>                     Bezier_polygon_vector ;
typedef Bezier_gps_traits::General_polygon_with_holes_2 Bezier_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Bezier_gps_traits>  Bezier_polygon_set ;

struct Bezier_helper
{
  static double get_approximated_endpoint_parameter( Bezier_point const& p, Bezier_curve const& curve, unsigned int xid  )
  {
    typedef Bezier_point::Originator_iterator Originator_iterator ;
    
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

  template<class OutputIterator>
  static void get_control_points ( Bezier_curve const& aCurve, bool aFwd, OutputIterator aOut )
  {
    int nc = aCurve.number_of_control_points() ;
    
    for ( int i = aFwd ? 0 : nc - 1 ; aFwd ? i < nc : i >= 0 ; aFwd ? ++ i : -- i )
    {
      *aOut++ = Linear_point( CGAL::to_double( aCurve.control_point(i).x() )
                            , CGAL::to_double( aCurve.control_point(i).y() )
                            ) ;
    }
  }  
  
  template<class OutputIterator>
  static void clip ( Bezier_curve const& aCurve, double aS, double aT, bool aFwd, OutputIterator aOut )
  {
    std::vector<Linear_point> lQ ;
    
    get_control_points(aCurve, aFwd, std::back_inserter(lQ) ) ;
    
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
  
  template<class OutputIterator>
  static void clip ( Bezier_X_monotone_curve const& aBXMC, OutputIterator aOut )
  {
    Bezier_curve const& lBC = aBXMC.supporting_curve();
  
    double lS = get_approximated_endpoint_parameter(aBXMC.source(), lBC, aBXMC.xid());
    double lT = get_approximated_endpoint_parameter(aBXMC.target(), lBC, aBXMC.xid());
    
    bool lFwd = lS <= lT ;
    
    double lMin = lFwd ? lS : 1.0 - lS ;
    double lMax = lFwd ? lT : 1.0 - lT ;
    
    clip(lBC, lMin, lMax, lFwd, aOut ); 
  }
} ;

struct Compute_bezier_X_monotone_cuve_bbox
{
  CGAL::Bbox_2 operator()( Bezier_X_monotone_curve const& aCurve ) const 
  {
    return aCurve.supporting_curve().bbox();
    
    std::vector<Linear_point> lQ ;
    
    Bezier_helper::get_control_points(aCurve.supporting_curve(), true, std::back_inserter(lQ) ) ;
    
    return CGAL::bbox_2(lQ.begin(),lQ.end());
  }
} ;

struct Draw_bezier_X_monotone_cuve
{
  template<class Path, class Converter>
  void operator()( Bezier_X_monotone_curve const& aBXMC, Path* aPath, Converter aConvert, int aIdx ) const 
  {
    std::vector<Linear_point> lQ ;
    Bezier_helper::clip(aBXMC,std::back_inserter(lQ));
    if ( lQ.size() == 2 )
    {
      if ( aIdx == 0 )
           aPath->moveTo( aConvert(lQ[0]) ) ;
      else aPath->lineTo( aConvert(lQ[0]) ) ; 
      
      aPath->lineTo( aConvert(lQ[1]) ) ;
    }
    else if ( lQ.size() == 3 )
    {
      if ( aIdx == 0 )
           aPath->moveTo ( aConvert(lQ[0]) ) ;
      else aPath->lineTo ( aConvert(lQ[0]) ) ; 
      
      aPath->quadTo( aConvert(lQ[1]), aConvert(lQ[2]) ) ; 
    }
    else if ( lQ.size() == 4 )
    {
      if ( aIdx == 0 )
           aPath->moveTo ( aConvert(lQ[0]) ) ;
      else aPath->lineTo ( aConvert(lQ[0]) ) ; 
      
      aPath->cubicTo( aConvert(lQ[1]), aConvert(lQ[2]), aConvert(lQ[3]) ) ; 
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
             aPath->moveTo(lP) ;
        else aPath->lineTo(lP) ; 
      }
    }
  }
} ;

typedef CGAL::Qt::GeneralPolygonSetGraphicsItem<Bezier_polygon_set,Compute_bezier_X_monotone_cuve_bbox,Draw_bezier_X_monotone_cuve> Bezier_GI;


#endif

#endif
