// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Laurent Saboret


#ifndef CGAL_ORIENT_NORMALS_WRT_CAMERAS_3_H
#define CGAL_ORIENT_NORMALS_WRT_CAMERAS_3_H

#include <CGAL/basic.h>
#include <CGAL/Oriented_normal_3.h>

#include <iterator>

CGAL_BEGIN_NAMESPACE


/// Orient a 3D point's normal w.r.t. the position of cameras
/// that reconstructed the point by photogrammetry.
template < class Gt, ///< Geometric traits class.
           class OrientedNormal_3,
           class InputIterator
>
void
orient_normal_wrt_cameras_3(const typename Gt::Point_3& p, ///< 3D point position
                            OrientedNormal_3& normal, ///< 3D point normal (in and out)
                            InputIterator first_camera,  ///< 3D point cameras
                            InputIterator beyond_camera) ///< 3D point cameras
{
    typedef typename Gt::FT       FT;
    typedef typename Gt::Point_3  Point;
    typedef OrientedNormal_3      Normal;
    typedef typename Gt::Vector_3 Vector;

    Vector n = normal.get_vector();

    //                                        ->
    //                                        cp     ->
    // Search for the camera position c / | ------ * n | is maximum
    //                                      ||cp||
    FT max_dot_product = 0;
    Point max_camera;
    for(InputIterator c = first_camera; c != beyond_camera; c++)
    {
    Vector cp = p - *c;
    FT norm_cp = std::sqrt(cp * cp);
    if (norm_cp != 0.0)
    {
      FT dot = (cp * n) / norm_cp;
      if (dot < 0)
        dot = -dot;

      if (max_dot_product < dot)
      {
         max_dot_product = dot;
         max_camera = *c;
      }
    }
    }

    //        ->             ->
    // Orient n backwards to cp
    if (max_dot_product > 0)
    {
      Vector cp = p - max_camera;
      FT dot = (cp * n);
      if (dot > 0)
         n = -n;
      normal = Normal(n, true /* oriented */);
    }
    else // if failure
    {
      normal = Normal(n, false /* non oriented */);
    }
    //// TEST: flag only inverted normals as oriented to see the result in 3D rendering
    //if (max_dot_product > 0)
    //{
    //  Vector cp = p - max_camera;
    //  FT dot = (cp * n);
    //  if (dot > 0)
    //    normal = Normal(-n, true /* oriented */);
    //  else
    //    normal = Normal(n, false /* non oriented */);
    //}
    //else // if failure
    //{
    //  normal = Normal(n, false /* non oriented */);
    //}
}


/// Orient the normals of the [first, beyond) range of vertices
/// w.r.t. the position of cameras
/// that reconstructed the points by photogrammetry.
///
/// Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - VertexCamerasMap is a model of boost::readable_property_map.

template < class VertexIterator, class VertexPointMap, class VertexNormalMap, class VertexCamerasMap >
void
orient_normals_wrt_cameras_3(VertexIterator first, ///< range of vertices
                             VertexIterator beyond, ///< range of vertices
                             VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
                             VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
                             VertexCamerasMap vertex_cameras_map) ///< property map VertexIterator -> pair of camera iterators
{
  typedef typename std::iterator_traits<VertexIterator>::value_type Vertex_type;
  typedef typename Vertex_type::Geom_traits Gt;

  // iterate over input points and orient normals
  for (VertexIterator it = first; it != beyond; it++)
  {
      orient_normal_wrt_cameras_3<Gt>(get(vertex_point_map,it),
                                      vertex_normal_map[it],
                                      get(vertex_cameras_map,it).first, get(vertex_cameras_map,it).second);
  }
}


CGAL_END_NAMESPACE

#endif // CGAL_ORIENT_NORMALS_WRT_CAMERAS_3_H

