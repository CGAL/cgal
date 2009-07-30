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

#ifndef CGAL_CIRCULAR_POLYGONAL_SAMPLER_2_H
#define CGAL_CIRCULAR_POLYGONAL_SAMPLER_2_H

#include <iterator>
#include <deque>

#include <CGAL/Coercion_traits.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace CGAL {

template<class CircularXMonotoneCurve, class OutputPointIterator>
void sample_circular_X_monotone_curve_2( CircularXMonotoneCurve const& aCXMC, OutputPointIterator aOut )
{
  typedef typename value_type_traits<OutputPointIterator>::type Output_point_type ;
  
  typedef std::vector<std::pair<double, double> > Samples_vector ;
  
  Samples_vector lSamples ;
  
  if(aCXMC->is_linear())
  {
    // when the curve is linear approximate will allways return
    // two pairs (for each endpoint) regardless the parameter of number
    // of points
    citr->approximate(std::back_inserter(lSamples), 0);
  }
  else
  {
    // circular arc
    double sx = CGAL::to_double(aCXMC->source().x());
    double tx = CGAL::to_double(aCXMC->target().x());
    int x_min;
    int x_max;
    if(aCXMC->is_directed_right())
    {
      x_min =  sx;
      x_max =  tx;
    }
    else
    {
      x_min = tx;
      x_max = sx;
    }
    
    const int n = ( x_max - x_min ) * 100 ;
    if (n > 0)
      aCXMC->approximate(std::back_inserter(lSamples), n);
  }
      
  for( typename Samples_vector::const_iterator it = lSamples.begin(); it != lSamples.end(); ++ it )
    *aOut ++ = Output_point_type(it->first, it->second);
}

template< class CircularPolygon, class Polygon>
void sample_circular_polygon_2( CircularPolygon const& aBP, Polygon& rPoly )
{
  for( typename CircularPolygon::Curve_const_iterator cit = aBP.curves_begin(); cit != aBP.curves_end(); ++ cit )
    sample_circular_X_monotone_curve_2(*cit, std::back_inserter(rPoly) );    
}

template<class PolygonWihHoles, class CircularPolygonWithHoles>
void sample_circular_polygon_with_holes_2( CircularPolygonWithHoles const& aCPWH, PolygonWihHoles& rPolyWithHoles )
{
  sample_circular_polygon_2(aCPWH.outer_boundary(), rPolyWithHoles.outer_boundary() ) ;
  
  for( typename CircularPolygonWithHoles::Hole_const_iterator hit = aCPWH.holes_begin(); hit != aCPWH.holes_end(); ++ hit )
  {
    typename PolygonWihHoles::Polygon_2 lHole ;
    sample_circular_polygon_2(*hit, lHole );    
    rPolyWithHoles.add_hole(lHole);
    
  }
}

template<class PolygonWithHoles>
struct Circular_polygon_with_holes_sampler
{
  typedef PolygonWithHoles result_type ;
  
  template<class CircularPolygonWithHoles>
  result_type operator()( CircularPolygonWithHoles const& aCPWH )
  {
    result_type rR ;
    sample_circular_polygon_with_holes_2(aCPWH, rR);
    return rR ;  
  }
  
} ;
} // namespace CGAL

#endif // CGAL_CIRCULAR_POLYGONAL_SAMPLER_2_H
