// Copyright (c) 2001,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_STATIC_FILTERS_H
#define CGAL_STATIC_FILTERS_H

// This kernel wrapper gathers optimized predicates written by hand, using
// a few steps of filtering.  It should work if the initial traits has
// cartesian coordinates which fit exactly in doubles.
//
// Purely static filters code has been removed, since it requires additional
// logic and is not plug'n play (requires users providing bounds).
// If it should be provided again, it should probably be separate.

#include <CGAL/basic.h>

#include <CGAL/Kernel/function_objects.h>
#include <CGAL/Cartesian/function_objects.h>

#include <CGAL/Static_filters/tools.h>
#include <CGAL/Static_filters/Orientation_2.h>
#include <CGAL/Static_filters/Orientation_3.h>
#include <CGAL/Static_filters/Side_of_oriented_circle_2.h>
#include <CGAL/Static_filters/Side_of_oriented_sphere_3.h>

// #include <CGAL/Static_filters/Coplanar_orientation_3.h>
// #include <CGAL/Static_filters/Coplanar_side_of_bounded_circle_3.h>

// TODO :
// - aim at obsoleting Filtered_exact, so that
//   Exact_predicates_inexact_constructions_kernel becomes Filtered_kernel.
// - Is Fixed_precision_nt now obsolete ?  If yes, deprecate it.
// - add more predicates :
//   - lexicographical comparisons
//   - left_turn (via generic adapter to orientation)
//   - power_tests
//   - others ?
// - benchmark on more algorithms.
// - improve fit_in_double() for other NTs (MP_Float, Lazy).
//   possibly promote it to NT requirement ?
// - try to automatize : have a struct a la Static_filter_error, with one part
//   which is runtime, and the other which can be constant-propagated by the
//   compiler.  g++ 4.0 should be able to cprop the second part...


CGAL_BEGIN_NAMESPACE

// The K_base argument is supposed to provide exact primitives.
template < typename K_base >
class Static_filters : public K_base
{
  typedef Static_filters<K_base>                    Self;

public:

  typedef CommonKernelFunctors::Left_turn_2<Self>   Left_turn_2;
  typedef CartesianKernelFunctors::Less_xy_2<Self>  Less_xy_2;
  typedef CartesianKernelFunctors::Less_yx_2<Self>  Less_yx_2;

  typedef SF_Orientation_2<K_base>                  Orientation_2;
  typedef SF_Orientation_3<K_base>                  Orientation_3;
  typedef SF_Side_of_oriented_circle_2<K_base>      Side_of_oriented_circle_2;
  typedef SF_Side_of_oriented_sphere_3<K_base>      Side_of_oriented_sphere_3;


  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(); }

  Less_xy_2
  less_xy_2_object() const
  { return Less_xy_2(); }

  Less_yx_2
  less_yx_2_object() const
  { return Less_yx_2(); }

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(); }

  Orientation_3
  orientation_3_object() const
  { return Orientation_3(); }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }

  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const
  { return Side_of_oriented_sphere_3(); }

  // The two following are for degenerate cases, so I'll update them later.
  //
  // typedef SF_Coplanar_orientation_3<Point_3, Orientation_2>
  //                                                   Coplanar_orientation_3;
  // typedef SF_Side_of_bounded_circle_3<Point_3>
  //                                         Coplanar_side_of_bounded_circle_3;

  // Coplanar_orientation_3
  // coplanar_orientation_3_object() const
  // { return Coplanar_orientation_3(); }

  // Coplanar_side_of_bounded_circle_3
  // coplanar_side_of_bounded_circle_3_object() const
  // { return Coplanar_side_of_bounded_circle_3(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_STATIC_FILTERS_H
