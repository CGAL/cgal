// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman     <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENV_PLANE_TRAITE_3_FUNCTIONS_H
#define CGAL_ENV_PLANE_TRAITE_3_FUNCTIONS_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/use.h>

namespace CGAL {

template <class K>
Object plane_half_plane_proj_intersection(const typename K::Plane_3 &h1, 
                                          const typename K::Plane_3 &h2,
                                          const typename K::Line_2  &l,
                                          const K& k)
{
  typedef typename K::Line_3     Line_3;
  typedef typename K::Line_2     Line_2;
  typedef typename K::Plane_3    Plane_3;

  // intersect the two planes
  Object h_obj = k.intersect_3_object()(h1, h2);
  if(h_obj.is_empty())
    return Object(); // no intersection at all (paralles planes)

  Plane_3 p;
  if(assign(p, h_obj))
    return Object();

  // if two planes are not parallel they must intersect at a 3D line
  Line_3 l3;
  CGAL_assertion_code(bool b =)
    assign(l3, h_obj);
  CGAL_assertion(b);

  const Line_2& proj_inter_line = project_xy(l3, k);

  return line_under_linear_constraint(proj_inter_line, l, k);
}

template <class K>
Object half_plane_half_plane_proj_intersection(const typename K::Plane_3 &h1,
                                               const typename K::Line_2  &l1,
                                               const typename K::Plane_3 &h2,
                                               const typename K::Line_2  &l2,
                                               const K& k)
{
  typedef typename K::Ray_2      Ray_2;
  typedef typename K::Line_2     Line_2;

  Object obj = plane_half_plane_proj_intersection(h1, h2, l2, k);
  if(obj.is_empty())
    return Object();

  Line_2 l;
  if(assign(l, obj))
    return line_under_linear_constraint(l, l1, k);

  Ray_2 ray;
  if(assign(ray, obj))
    return ray_under_linear_constraint(ray, l1, k);

  CGAL_error(); // doesnt suppose to reach here
  return Object();
}

template <class K>
typename K::Line_2 project_xy(const typename K::Line_3& l,
                              const K& k)
{
  typedef typename K::Vector_3   Vector_3;
  typedef typename K::Vector_2   Vector_2;
  typedef typename K::Point_3    Point_3;
  typedef typename K::Point_2    Point_2;

  Vector_3 vec3 = k.construct_vector_3_object()(l);
  Vector_2 vec2(vec3.x(), vec3.y());
  Point_3 p3 = k.construct_point_on_3_object()(l, 0);
  Point_2 p2(p3.x(), p3.y());
  return k.construct_line_2_object()(p2, vec2);
}


// l1 is a line, l2 is a linear constraint (determined by the direction
// of the line).
template <class K>
Object line_under_linear_constraint(const typename K::Line_2& l1,
                                    const typename K::Line_2& l2,
                                    const K& k)
{
  typedef typename K::Ray_2         Ray_2;
  typedef typename K::Line_2        Line_2;
  typedef typename K::Vector_2      Vector_2;
  typedef typename K::Point_2       Point_2;

  Object obj = k.intersect_2_object()(l1, l2);
  Point_2 p;
  if(assign(p, obj))
  {
    const Vector_2& vec = k.construct_vector_2_object()(l1);
    const Point_2& s = k.construct_translated_point_2_object()(p, vec);
    const Ray_2& ray = k.construct_ray_2_object()(p, s);
    Oriented_side side = k.oriented_side_2_object()(l2, s);
    if(side == ON_NEGATIVE_SIDE)
    {
      return make_object(k.construct_opposite_ray_2_object()(ray));
    }

    CGAL_assertion(side == ON_POSITIVE_SIDE); //the two lines are not parallel
    return make_object(ray);
  }

  if(obj.is_empty()) // the two lines are parallel
  {
    const Point_2& s = k.construct_point_on_2_object()(l1, 0);
    Oriented_side side = k.oriented_side_2_object()(l2, s);

    if(side == ON_NEGATIVE_SIDE)
      return Object();

    CGAL_assertion(side == ON_POSITIVE_SIDE); // the two lines are parallel
    return make_object(l1);   
  }
 
  // the two lines overlap
  CGAL_USE_TYPE(Line_2);
  CGAL_assertion_code(Line_2 dummy;);
  CGAL_assertion_code(bool b =  assign(dummy, obj););
  CGAL_assertion(b);
  
  return make_object(l1);
}

template <class K>
Object ray_under_linear_constraint(const typename K::Ray_2&  ray,
                                   const typename K::Line_2& l,
                                   const K& k)
{
  typedef typename K::Vector_2      Vector_2;
  typedef typename K::Point_2       Point_2;

  const Point_2& s = k.construct_point_on_2_object()(ray, 0);
  Oriented_side side = k.oriented_side_2_object()(l, s);
  Object obj = k.intersect_2_object()(ray, l);
  if(obj.is_empty())
  {
    if(side == ON_NEGATIVE_SIDE)
      return Object();

    CGAL_assertion(side == ON_POSITIVE_SIDE);
    return make_object(ray);
  }

  Point_2 p;
  if(assign(p, obj))
  {
    if(side == ON_POSITIVE_SIDE)
      return make_object(k.construct_segment_2_object()(s, p));

    Vector_2 vec = k.construct_vector_2_object()(ray);
    if(side == ON_NEGATIVE_SIDE)
      return make_object(k.construct_ray_2_object()(p, vec));

    CGAL_assertion(side == ON_ORIENTED_BOUNDARY);
    const Vector_2& vec_of_l = k.construct_vector_2_object()(l);

    Orientation orient = k.orientation_2_object()(vec_of_l, vec);
    if(orient == LEFT_TURN)
      return make_object(ray);

    CGAL_assertion(orient == RIGHT_TURN);
    return make_object(s);
  }

  // the ray and the line overlap
  return make_object(ray);
}
 

} //namespace CGAL

#endif
