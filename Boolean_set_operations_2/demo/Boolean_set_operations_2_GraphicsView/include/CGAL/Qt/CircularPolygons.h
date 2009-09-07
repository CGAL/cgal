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
  template<class Circular_X_monotone_curve, class Path, class Converter>
  void operator()( Circular_X_monotone_curve const& curve, Path& aPath, Converter aConvert, int aIdx ) const 
  {
    typedef Point_2< Simple_cartesian<double> > Linear_point ;
    
    Linear_point lS( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
    Linear_point lT( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
    
    if ( aIdx == 0 )
         aPath.moveTo( aConvert( lS ) ) ;
    else aPath.lineTo( aConvert( lS ) ) ;
    
    aPath.lineTo( aConvert( lT ) ) ;
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
