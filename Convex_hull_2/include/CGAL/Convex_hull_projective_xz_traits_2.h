// Copyright (c) 2001  Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Susan Hert

#ifndef CGAL_CONVEX_HULL_PROJECTIVE_XZ_TRAITS_2_H
#define CGAL_CONVEX_HULL_PROJECTIVE_XZ_TRAITS_2_H

#include <CGAL/license/Convex_hull_2.h>


#define CGAL_DEPRECATED_HEADER "<CGAL/Convex_hull_projective_xz_traits_2.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Projection_traits_xz_3.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/function_objects.h>

namespace CGAL {

template <class Point_3>
class Less_xy_plane_xz_2 
{
public:
   typedef bool           result_type;

   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.x(), p.z(), q.x(), q.z()) == SMALLER;
   }
};

template <class Point_3>
class Equal_xy_plane_xz_2 
{
public:
   typedef bool           result_type;

   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.x(), p.z(), q.x(), q.z()) == EQUAL;
   }
};

template <class Point_3>
class Less_yx_plane_xz_2 
{
public:
   typedef bool           result_type;

   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.z(), p.x(), q.z(), q.x()) == SMALLER;
   }
};

template <class Point_3>
class Left_turn_plane_xz_2 
{
public:
   typedef bool           result_type;

   bool 
   operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
   { 
    return orientationC2(p.x(), p.z(), q.x(), q.z(), r.x(), r.z()) 
                                                          == LEFT_TURN;
   }
};


template <class Point_3>
class Less_dist_to_line_plane_xz_2
{
public:
   typedef bool           result_type;

   bool
   operator()(const Point_3& p, const Point_3& q,
              const Point_3& r, const Point_3& s) const
   {
      Comparison_result
         res = cmp_signed_dist_to_lineC2(p.x(), p.z(), q.x(), q.z(),
                                         r.x(), r.z(), s.x(), s.z());
      if ( res == LARGER )
         return false;
      else if ( res == SMALLER )
         return true;
      else
         return compare_lexicographically_xyC2(r.x(), r.z(), s.x(), s.z())
             == SMALLER;
   }
};


template <class Point_3>
class Less_rotate_ccw_plane_xz_2
{
public:

   typedef bool         result_type;

   bool
   operator()(const Point_3& r, const Point_3& p, const Point_3& q) const
   {
      Orientation orient =
               orientationC2(r.x(), r.z(), p.x(), p.z(), q.x(), q.z());
      if ( orient ==  LEFT_TURN )
         return true;
      else if ( orient == RIGHT_TURN )
         return false;
      else
      {
         if (p.x() == r.x() && p.z() == r.z()) return false;
         if (q.x() == r.x() && q.z() == r.z()) return true;
         if (p.x() == q.x() && p.z() == q.z()) return false;
         return
            collinear_are_ordered_along_lineC2(r.x(), r.z(),
                                               q.x(), q.z(), p.x(), p.z());
      }
   }

};



template <class Point_3>
class Convex_hull_projective_xz_traits_2 
{
public:
    typedef Point_3                             Point_2;
    typedef Less_xy_plane_xz_2<Point_3>         Less_xy_2;
    typedef Equal_xy_plane_xz_2<Point_3>        Equal_2;    
    typedef Less_yx_plane_xz_2<Point_3>         Less_yx_2;
    typedef Left_turn_plane_xz_2<Point_3>       Left_turn_2;
    typedef Less_rotate_ccw_plane_xz_2<Point_3> Less_rotate_ccw_2;
    typedef Less_dist_to_line_plane_xz_2<Point_3> 
                                                Less_signed_distance_to_line_2;
    Less_xy_2
    less_xy_2_object() const
    {  return Less_xy_2(); }
    
    Equal_2
    equal_2_object() const
    {  return Equal_2(); }    

    Less_yx_2
    less_yx_2_object() const
    {  return Less_yx_2(); }

    Left_turn_2
    left_turn_2_object() const
    {  return Left_turn_2(); }

    Less_rotate_ccw_2
    less_rotate_ccw_2_object() const
    {  return Less_rotate_ccw_2(); }

    Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object() const
    {  return Less_signed_distance_to_line_2(); }

};

}

#endif // CGAL_CONVEX_HULL_PROJECTIVE_XZ_TRAITS_2_H
