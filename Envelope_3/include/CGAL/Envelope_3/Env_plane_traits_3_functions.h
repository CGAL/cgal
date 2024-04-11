// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman     <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENV_PLANE_TRAITE_3_FUNCTIONS_H
#define CGAL_ENV_PLANE_TRAITE_3_FUNCTIONS_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/basic.h>
#include <CGAL/use.h>

namespace CGAL {

template <class K>
std::optional< std::variant<typename K::Line_2, typename K::Ray_2, typename K::Segment_2, typename K::Point_2> >
plane_half_plane_proj_intersection(const typename K::Plane_3 &h1,
                                   const typename K::Plane_3 &h2,
                                   const typename K::Line_2  &l,
                                   const K& k)
{
  typedef typename K::Line_3     Line_3;
  typedef typename K::Line_2     Line_2;
  typedef typename K::Plane_3    Plane_3;

  // intersect the two planes
  auto h_obj = k.intersect_3_object()(h1, h2);
  if(h_obj == std::nullopt)
    return std::nullopt; // no intersection at all (parallel planes)

  Plane_3 p;
  if(std::get_if<Plane_3>(&(*h_obj))==nullptr)
    return std::nullopt;

  // if two planes are not parallel they must intersect at a 3D line
  const Line_3* l3 = std::get_if<Line_3>(&(*h_obj));
  CGAL_assertion(l3!=nullptr);

  const Line_2& proj_inter_line = project_xy(*l3, k);

  return line_under_linear_constraint(proj_inter_line, l, k); //LR
}

template <class K>
std::optional< std::variant<typename K::Line_2, typename K::Ray_2, typename K::Segment_2, typename K::Point_2> >
half_plane_half_plane_proj_intersection(const typename K::Plane_3 &h1,
                                        const typename K::Line_2  &l1,
                                        const typename K::Plane_3 &h2,
                                        const typename K::Line_2  &l2,
                                        const K& k)
{
  typedef typename K::Ray_2      Ray_2;
  typedef typename K::Line_2     Line_2;

  auto obj = plane_half_plane_proj_intersection(h1, h2, l2, k);
  if(obj == std::nullopt)
    return std::nullopt;

  if(const Line_2* line = std::get_if<Line_2>(&(*obj)))
    return line_under_linear_constraint(*line, l1, k);

  if(const Ray_2* ray = std::get_if<Ray_2>(&(*obj)))
    return ray_under_linear_constraint(*ray, l1, k);

  CGAL_error(); // doesn't suppose to reach here
  return std::nullopt;
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
std::optional< std::variant<typename K::Line_2, typename K::Ray_2, typename K::Segment_2, typename K::Point_2> >
line_under_linear_constraint(const typename K::Line_2& l1,
                             const typename K::Line_2& l2,
                             const K& k)
{
  typedef typename K::Ray_2         Ray_2;
  typedef typename K::Vector_2      Vector_2;
  typedef typename K::Point_2       Point_2;

  auto obj = k.intersect_2_object()(l1, l2);

  if(obj == std::nullopt)// the two lines are parallel
  {
    const Point_2& s = k.construct_point_on_2_object()(l1, 0);
    Oriented_side side = k.oriented_side_2_object()(l2, s);

    if(side == ON_NEGATIVE_SIDE)
      return std::nullopt;

    CGAL_assertion(side == ON_POSITIVE_SIDE); // the two lines are parallel
    return l1;
  }

  if(const Point_2* p = std::get_if<Point_2>(&(*obj)))
  {
    const Vector_2& vec = k.construct_vector_2_object()(l1);
    const Point_2& s = k.construct_translated_point_2_object()(*p, vec);
    const Ray_2& ray = k.construct_ray_2_object()(*p, s);
    Oriented_side side = k.oriented_side_2_object()(l2, s);
    if(side == ON_NEGATIVE_SIDE)
    {
      return k.construct_opposite_ray_2_object()(ray);
    }

    CGAL_assertion(side == ON_POSITIVE_SIDE); //the two lines are not parallel
    return ray;
  }

  // the two lines overlap
  return l1;
}

template <class K>
std::optional< std::variant<typename K::Line_2, typename K::Ray_2, typename K::Segment_2, typename K::Point_2> >
ray_under_linear_constraint(const typename K::Ray_2&  ray,
                            const typename K::Line_2& l,
                            const K& k)
{
  typedef typename K::Vector_2      Vector_2;
  typedef typename K::Point_2       Point_2;

  const Point_2& s = k.construct_point_on_2_object()(ray, 0);
  Oriented_side side = k.oriented_side_2_object()(l, s);
  auto obj = k.intersect_2_object()(ray, l);
  if(obj == std::nullopt)
  {
    if(side == ON_NEGATIVE_SIDE)
      return std::nullopt;

    CGAL_assertion(side == ON_POSITIVE_SIDE);
    return ray;
  }

  if(const Point_2* p = std::get_if<Point_2>(&(*obj)))
  {
    if(side == ON_POSITIVE_SIDE)
      return k.construct_segment_2_object()(s, *p);

    Vector_2 vec = k.construct_vector_2_object()(ray);
    if(side == ON_NEGATIVE_SIDE)
      return k.construct_ray_2_object()(*p, vec);

    CGAL_assertion(side == ON_ORIENTED_BOUNDARY);
    const Vector_2& vec_of_l = k.construct_vector_2_object()(l);

    Orientation orient = k.orientation_2_object()(vec_of_l, vec);
    if(orient == LEFT_TURN)
      return ray;

    CGAL_assertion(orient == RIGHT_TURN);
    return s;
  }

  // the ray and the line overlap
  return ray;
}


} //namespace CGAL

#endif
