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

#ifndef CGAL_QT_CIRCULAR_POLYGONS_H
#define CGAL_QT_CIRCULAR_POLYGONS_H

#include <CGAL/Qt/PiecewiseSetGraphicsItem.h>

namespace CGAL {

namespace Qt {

struct Circular_X_monotone_bbox
{
  template<class X_monotone_circle_segment_2>
  CGAL::Bbox_2 operator()( X_monotone_circle_segment_2 const& aC ) const 
  {
    return aC.bbox();
  }
} ;

struct Circular_bbox
{
  template<class Circle_segment_2>
  CGAL::Bbox_2 operator()( Circle_segment_2 const& aC ) const 
  {
    double x_min = to_double(aC.source().x());
    double x_max = to_double(aC.target().x());   
    double y_min = to_double(aC.source().y()); 
    double y_max = to_double(aC.target().y());
    if(x_min > x_max)
    {
      std::swap(x_min, x_max);
      std::swap(y_min, y_max);
    }

    if(y_min > y_max)
      std::swap(y_min, y_max);

    if(aC.is_circular())
    {
      typedef typename Circle_segment_2::Circle_2 Circle_2 ;

      const Circle_2& circ = aC.supporting_circle();

      return circ.bbox();
    }
    
    return Bbox_2(x_min, y_min, x_max, y_max);
  }
} ;

struct Draw_circular_X_monotone_curve
{
  template<class X_monotone_circle_segment_2, class Path>
  void operator()( X_monotone_circle_segment_2 const& curve, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
    
    typedef Point_2<Linear_kernel> Linear_point ;
    
    typedef Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    if ( curve.is_circular() )
    {
      typename X_monotone_circle_segment_2::Circle_2        const& circ   = curve.supporting_circle();
      typename X_monotone_circle_segment_2::Kernel::Point_2 const& center = circ.center();
      typename X_monotone_circle_segment_2::Point_2         const& source = curve.source();
      typename X_monotone_circle_segment_2::Point_2         const& target = curve.target();

      double sx = to_double(source.x());
      double sy = to_double(source.y());
      double tx = to_double(target.x());
      double ty = to_double(target.y());
      double cx = to_double(center.x());
      double cy = to_double(center.y());

      bool degenerate = ( sx == tx ) && ( sy == ty ) ;

      if ( !degenerate )
      {
        double sdy = sy - cy ;
        double sdx = sx - cx ;
        double tdy = ty - cy ;
        double tdx = tx - cx ;

        double asource = std::atan2(sdy, sdx ); 
        double atarget = std::atan2(tdy, tdx );

        if( asource < 0.0)
          asource += 2 * CGAL_PI;

        if( atarget <= 0.0)
          atarget += 2 * CGAL_PI;

        if ( atarget  < asource )
          atarget += 2 * CGAL_PI;

        double aspan = atarget - asource ;  

        const double to_deg = 180.0/CGAL_PI;

        Orientation lO = curve.orientation() ;

        if ( lO == CLOCKWISE )
          aspan = 2.0 * CGAL_PI - aspan ;

        if ( aIdx == 0 ) 
             aPath.moveTo(sx,sy) ;
        else aPath.lineTo(sx,sy) ;

        QRectF bbox = convert(circ.bbox()) ;

        double dasource = std::atan2(-sdy, sdx ) * to_deg ;

        double daspan  = aspan * to_deg * ( lO == COUNTERCLOCKWISE ? -1.0 : +1.0) ;

        // This is to prevent the approximations to turn a tiny arc into a full circle by
        // inverting the relative ordering of the start, target angles.
        // We use the fact that an X-monotone arc can never span an angle greater than PI.
        if ( daspan < 270 )
        {
          aPath.arcTo(bbox , dasource, daspan );    
        }
        else
        {
         if ( aIdx == 0 ) 
            aPath.moveTo(sx,sy) ;
          aPath.lineTo(sx,sy) ;
        }
      }
      else
      {
        if ( aIdx == 0 ) 
          aPath.moveTo(sx,sy) ;
        aPath.lineTo(sx,sy) ;
      }
    }
    else
    {
      Linear_point lS( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
      Linear_point lT( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
      
      if ( aIdx == 0 ) 
           aPath.moveTo( convert( lS ) ) ;
      else aPath.lineTo( convert( lS ) ) ;

      aPath.lineTo( convert( lT ) ) ;
    }
  }
} ;


struct Draw_circular_curve
{
  template<class Circle_segment_2, class Path>
  void operator()( Circle_segment_2 const& curve, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
    
    typedef Point_2<Linear_kernel> Linear_point ;
    
    typedef Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    if ( curve.is_circular() )
    {
      typename Circle_segment_2::Circle_2        const& circ   = curve.supporting_circle();
      typename Circle_segment_2::Kernel::Point_2 const& center = circ.center();
      typename Circle_segment_2::Point_2         const& source = curve.source();
      typename Circle_segment_2::Point_2         const& target = curve.target();

      double sx = to_double(source.x());
      double sy = to_double(source.y());
      double tx = to_double(target.x());
      double ty = to_double(target.y());
      double cx = to_double(center.x());
      double cy = to_double(center.y());

      bool degenerate = ( sx == tx ) && ( sy == ty ) ;

      if ( !degenerate )
      {
        double sdy = sy - cy ;
        double sdx = sx - cx ;
        double tdy = ty - cy ;
        double tdx = tx - cx ;

        double asource = std::atan2(sdy, sdx ); 
        double atarget = std::atan2(tdy, tdx );

        if( asource < 0.0)
          asource += 2 * CGAL_PI;

        if( atarget <= 0.0)
          atarget += 2 * CGAL_PI;

        if ( atarget  < asource )
          atarget += 2 * CGAL_PI;

        double aspan = atarget - asource ;  

        const double to_deg = 180.0/CGAL_PI;

        Orientation lO = curve.orientation() ;

        if ( lO == CLOCKWISE )
          aspan = 2.0 * CGAL_PI - aspan ;
 
        if ( aIdx == 0 ) 
             aPath.moveTo(sx,sy) ;
        else aPath.lineTo(sx,sy) ;

        QRectF bbox = convert(circ.bbox()) ;

        double dasource = std::atan2(-sdy, sdx ) * to_deg ;

        double daspan  = aspan * to_deg * ( lO == COUNTERCLOCKWISE ? -1.0 : +1.0) ;

        aPath.arcTo(bbox , dasource, daspan );    
      }
    }
    else
    {
      Linear_point lS( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
      Linear_point lT( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
      
      if ( aIdx == 0 ) 
           aPath.moveTo( convert( lS ) ) ;
      else aPath.lineTo( convert( lS ) ) ;

      aPath.lineTo( convert( lT ) ) ;
    }
  }
} ;

template<class Circular_boundary_pieces>
class Circular_boundary_pieces_graphics_item : public Boundary_pieces_graphics_item<Circular_boundary_pieces,Draw_circular_curve,Circular_bbox>
{
  typedef Boundary_pieces_graphics_item<Circular_boundary_pieces,Draw_circular_curve,Circular_bbox> Base ;
  
public :

  Circular_boundary_pieces_graphics_item( Circular_boundary_pieces* aPieces ) : Base(aPieces) {}
} ;
    
template<class Circular_boundary>
class Circular_boundary_graphics_item : public Piecewise_boundary_graphics_item<Circular_boundary,Draw_circular_X_monotone_curve,Circular_X_monotone_bbox>
{
  typedef Piecewise_boundary_graphics_item<Circular_boundary,Draw_circular_X_monotone_curve,Circular_X_monotone_bbox> Base ;
  
public :

  Circular_boundary_graphics_item( Circular_boundary* aBoundary ) : Base(aBoundary) {}
} ;

template<class Circular_region>
class Circular_region_graphics_item : public Piecewise_region_graphics_item<Circular_region,Draw_circular_X_monotone_curve,Circular_X_monotone_bbox>
{

  typedef Piecewise_region_graphics_item<Circular_region,Draw_circular_X_monotone_curve,Circular_X_monotone_bbox> Base ;
  
public:

  Circular_region_graphics_item(Circular_region* aRegion ) : Base(aRegion) {}  
} ;

template<class Circular_set>
class Circular_set_graphics_item : public Piecewise_set_graphics_item<Circular_set,Draw_circular_X_monotone_curve,Circular_X_monotone_bbox>
{

  typedef Piecewise_set_graphics_item<Circular_set,Draw_circular_X_monotone_curve,Circular_X_monotone_bbox> Base ;
  
public:

  Circular_set_graphics_item(Circular_set* aSet) : Base(aSet) {}
} ;


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CIRCULAR_POLYGONS_H
