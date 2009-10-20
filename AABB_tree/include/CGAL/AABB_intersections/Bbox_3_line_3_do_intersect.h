// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
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
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez


#ifndef CGAL_LINE_3_BBOX_3_DO_INTERSECT_H
#define CGAL_LINE_3_BBOX_3_DO_INTERSECT_H

#include <CGAL/Line_3.h>
#include <CGAL/Bbox_3.h>

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

CGAL_BEGIN_NAMESPACE

namespace internal {

  template <class K>
  bool do_intersect(const typename K::Line_3& line, 
    const CGAL::Bbox_3& bbox,
    const K&)
  {
    typedef typename K::FT FT;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;

    const Point source = line.point(0);

    Point parameters[2];
    parameters[0] = Point(bbox.xmin(), bbox.ymin(), bbox.zmin());
    parameters[1] = Point(bbox.xmax(), bbox.ymax(), bbox.zmax());

    const Vector direction = line.to_vector();
    
    // We don't care about values
    FT tmin(0.0);
    FT tmax(0.0);
    
    bool is_dx_null = false;
    bool is_dy_null = false;
    
    if ( ! CGAL_NTS is_zero(direction.x()) )
    {
      FT inv_direction_x = FT(1)/direction.x();
      const int sign_x = inv_direction_x < FT(0);
      
      tmin = (parameters[sign_x].x() - source.x()) * inv_direction_x;
      tmax = (parameters[1-sign_x].x() - source.x()) * inv_direction_x;
    }
    else
    {
      // premature exit if x value of line is outside bbox
      if ( source.x() < parameters[0].x() || source.x() > parameters[1].x() )
        return false;
      
      is_dx_null = true;
    }
    
    if ( ! CGAL_NTS is_zero(direction.y()) )
    {
      FT inv_direction_y = FT(1)/direction.y();
      const int sign_y = inv_direction_y < FT(0);
      
      const FT tymin = (parameters[sign_y].y() - source.y()) * inv_direction_y;
      const FT tymax = (parameters[1-sign_y].y() - source.y()) * inv_direction_y;
      
      if ( !is_dx_null )
      {
        if(tmin > tymax || tymin > tmax) 
          return false;
        
        if(tymin > tmin)
          tmin = tymin;
        
        if(tymax < tmax)
          tmax = tymax;
      }
      else
      {
        tmin = tymin;
        tmax = tymax;
      }
    }
    else
    {
      // premature exit if y value of line is outside bbox
      if ( source.y() < parameters[0].y() || source.y() > parameters[1].y() )
        return false;
      
      is_dy_null = true;
    }
    
    if ( ! CGAL_NTS is_zero(direction.z()) )
    {
      FT inv_direction_z = FT(1)/direction.z();
      const int sign_z = inv_direction_z < FT(0);
      
      FT tzmin = (parameters[sign_z].z() - source.z()) * inv_direction_z;
      FT tzmax = (parameters[1-sign_z].z() - source.z()) * inv_direction_z;
      
      if( (!is_dx_null || !is_dy_null) && (tmin > tzmax || tzmin > tmax) ) 
        return false;
    }
    else
    {
      // premature exit if z value of line is outside bbox
      if ( source.z() < parameters[0].z() || source.z() > parameters[1].z() )
        return false;
    }
    
    return true;
  }

} // namespace internal

template <class K>
bool do_intersect(const CGAL::Line_3<K>& line, 
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(line, bbox);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox, 
		  const CGAL::Line_3<K>& line)
{
  return typename K::Do_intersect_3()(line, bbox);
}

CGAL_END_NAMESPACE

#endif  // CGAL_LINE_3_BBOX_3_DO_INTERSECT_H


