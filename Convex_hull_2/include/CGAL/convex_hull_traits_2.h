// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Stefan Schirra

// This file's name must begin with a lower-case letter for backward 
// compatability.  Unfortunately, you can't have a file that differs only 
// in capitalization on the Windows platforms.

#ifndef CGAL_CONVEX_HULL_TRAITS_2_H
#define CGAL_CONVEX_HULL_TRAITS_2_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/distance_predicates_2.h>

namespace CGAL {

template <class K_>
class Convex_hull_traits_2 : public K_
{
public:
  typedef K_                                 K;
  typedef typename K::Point_2                Point_2;    
  typedef typename K::Less_xy_2              Less_xy_2;
  typedef typename K::Less_yx_2              Less_yx_2;
  typedef typename K::Less_signed_distance_to_line_2  
                                         Less_signed_distance_to_line_2;
  typedef typename K::Less_rotate_ccw_2      Less_rotate_ccw_2;
  typedef typename K::Left_turn_2             Left_turn_2;
  typedef typename K::Equal_2                Equal_2;
  typedef typename K::Segment_2              Segment_2;    
  
  Less_xy_2
  less_xy_2_object() const 
  { return Less_xy_2(); } 

  Less_yx_2
  less_yx_2_object() const 
  { return Less_yx_2(); } 

  Less_signed_distance_to_line_2
  less_signed_distance_to_line_2_object( ) const
  { return Less_signed_distance_to_line_2( ); } 

  Less_rotate_ccw_2
  less_rotate_ccw_2_object() const
  { return Less_rotate_ccw_2(); }

  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(); }

  Equal_2
  equal_2_object() const
  { return Equal_2(); }
};


template <class K>
class convex_hull_traits_2 : public Convex_hull_traits_2<K>
{};

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_2_H
