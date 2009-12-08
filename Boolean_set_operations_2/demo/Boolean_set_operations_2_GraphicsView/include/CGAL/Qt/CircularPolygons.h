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

#ifndef CGAL_QT_CIRCULAR_POLYGONS_H
#define CGAL_QT_CIRCULAR_POLYGONS_H

#include <CGAL/Qt/PiecewiseSetGraphicsItem.h>

namespace CGAL {

namespace Qt {

struct Circular_X_monotone_bbox
{
  template<class Circular_X_monotone_curve>
  CGAL::Bbox_2 operator()( Circular_X_monotone_curve const& aCXMC ) const 
  {
    return aCXMC.bbox();
  }
} ;

struct Draw_circular_X_monotone_curve
{
  template<class Circular_X_monotone_curve, class Path>
  void operator()( Circular_X_monotone_curve const& curve, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
    
    typedef Point_2<Linear_kernel> Linear_point ;
    
    typedef Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    if ( curve.is_circular() )
    {
      typename Circular_X_monotone_curve::Circle_2        const& circ   = curve.supporting_circle();
      typename Circular_X_monotone_curve::Kernel::Point_2 const& center = circ.center();
      typename Circular_X_monotone_curve::Point_2         const& source = curve.source();
      typename Circular_X_monotone_curve::Point_2         const& target = curve.target();

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

        // This is to prevent the approximations to turn a tiny arc into a full circle by
        // inverting the relative ordering of the start, target angles.
        // We use the fact that an X-monotone arc can never span an angle greater than PI.
        if ( aspan < 3.0 * CGAL_PI / 2.0 )
        {
          const double to_deg = 180.0/CGAL_PI;

          Orientation lO = curve.orientation() ;

          if ( aIdx == 0 ) 
               aPath.moveTo(sx,sy) ;
          else aPath.lineTo(sx,sy) ;

          QRectF bbox = convert(circ.bbox()) ;

          double dasource = std::atan2(-sdy, sdx ) * to_deg ;

          double daspan  = aspan * to_deg * ( lO == COUNTERCLOCKWISE ? -1.0 : +1.0) ;

          aPath.arcTo(bbox , dasource, daspan );    
        }
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
