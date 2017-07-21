// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_INTERNAL_STATIC_FILTERS_H
#define CGAL_INTERNAL_STATIC_FILTERS_H

// This kernel wrapper gathers optimized predicates written by hand, using
// a few steps of filtering.  It should work if the initial traits has
// cartesian coordinates which fit exactly in doubles.
//
// Purely static filters code has been removed, since it requires additional
// logic and is not plug'n play (requires users providing bounds).
// If it should be provided again, it should probably be separate.

#include <CGAL/config.h>

#include <CGAL/Kernel/function_objects.h>
#include <CGAL/Cartesian/function_objects.h>

#include <CGAL/internal/Static_filters/tools.h>
#include <CGAL/internal/Static_filters/Orientation_2.h>
#include <CGAL/internal/Static_filters/Orientation_3.h>
#include <CGAL/internal/Static_filters/Collinear_3.h>

// for static filters added nov./dec. 2011
#ifdef CGAL_DISABLE_STATIC_FILTERS_ADDED_2011
#  define CGAL_NO_EQUAL_3_STATIC_FILTERS 1
#  define CGAL_NO_COMPARE_X_2_STATIC_FILTERS 1
#  define CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS 1
#  define CGAL_NO_ANGLE_3_STATIC_FILTERS 1
#  define CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS 1
#endif // CGAL_DISABLE_STATIC_FILTERS_ADDED_2011


#ifndef CGAL_NO_EQUAL_3_STATIC_FILTERS
#  include <CGAL/internal/Static_filters/Equal_3.h>
#  include <CGAL/internal/Static_filters/Equal_2.h>
#endif // NOT CGAL_NO_EQUAL_3_STATIC_FILTERS

#ifndef CGAL_NO_COMPARE_X_2_STATIC_FILTERS
#  include <CGAL/internal/Static_filters/Compare_x_2.h>
#  include <CGAL/internal/Static_filters/Compare_y_2.h>
#endif // NOT CGAL_NO_COMPARE_X_2_STATIC_FILTERS

#ifndef CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS
#  include <CGAL/internal/Static_filters/Is_degenerate_3.h>
#endif // NOT CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS

#ifndef CGAL_NO_ANGLE_3_STATIC_FILTERS
#  include <CGAL/internal/Static_filters/Angle_3.h>
#endif // NOT CGAL_NO_ANGLE_3_STATIC_FILTERS

#ifndef CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS
#  include <CGAL/internal/Static_filters/Do_intersect_3.h>
#endif // NOT NOT CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS

#include <CGAL/internal/Static_filters/Compare_y_at_x_2.h>
#include <CGAL/internal/Static_filters/Side_of_oriented_circle_2.h>
#include <CGAL/internal/Static_filters/Side_of_oriented_sphere_3.h>
#include <CGAL/internal/Static_filters/Compare_squared_radius_3.h>
#include <CGAL/internal/Static_filters/Compare_weighted_squared_radius_3.h>
#include <CGAL/internal/Static_filters/Power_side_of_oriented_power_sphere_3.h>


// #include <CGAL/internal/Static_filters/Coplanar_orientation_3.h>
// #include <CGAL/internal/Static_filters/Coplanar_side_of_bounded_circle_3.h>

// TODO :
// - add more predicates :
//   - all lexicographical comparisons
//   - left_turn (via generic adapter to orientation)
//   - others ?
// - benchmark on more algorithms.
// - improve fit_in_double() for other NTs (MP_Float, Lazy). cf tools.h.
// - try to automatize : have a struct a la Static_filter_error, with one part
//   which is runtime, and the other which can be constant-propagated by the
//   compiler.  g++ 4.0 should be able to cprop the second part...


// Note about the second parameter of Static_filters<K,bool>:
// - if the access to Cartesian exact coordinates is cheap
//   (Simple_cartesian, Cartesian), then one can implement predicates that
//   just compare coordinates without filtering, using unfiltered
//   predicates defined in the namespace CartesianKernelFunctors.
// 
// - in the case of Lazy_kernel, where the access to p.x(), for a point p,
//   triggers the construction of a Lazy_exact_nt object, one does not want
//   to use the functors from the namespace CartesianKernelFunctors.

namespace CGAL { namespace internal {

// Here is the case when has_cheap_access_to_cartesian_coordinates is
// false, used by Lazy_kernel
// The K_base argument is supposed to provide exact primitives.
template < typename K_base, 
           bool has_cheap_access_to_cartesian_coordinates = true>
class Static_filters : public K_base {


  typedef Static_filters<K_base, 
                         has_cheap_access_to_cartesian_coordinates>         Self;

public:
#ifndef CGAL_NO_EQUAL_3_STATIC_FILTERS
  typedef Static_filters_predicates::Equal_2<K_base>                        Equal_2;
  typedef Static_filters_predicates::Equal_3<K_base>                        Equal_3;
#endif // NOT CGAL_NO_EQUAL_3_STATIC_FILTERS

#ifndef CGAL_NO_COMPARE_X_2_STATIC_FILTERS
  typedef Static_filters_predicates::Compare_x_2<K_base>                    Compare_x_2;
  typedef Static_filters_predicates::Compare_y_2<K_base>                    Compare_y_2;
#endif // NOT CGAL_NO_COMPARE_X_2_STATIC_FILTERS

#ifndef CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS
  typedef Static_filters_predicates::Is_degenerate_3<K_base, Self>          Is_degenerate_3;
#endif // NOT CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS
  typedef Static_filters_predicates::Orientation_2<K_base>                  Orientation_2;
  typedef Static_filters_predicates::Orientation_3<K_base>                  Orientation_3;
#ifndef CGAL_NO_ANGLE_3_STATIC_FILTERS

  typedef Static_filters_predicates::Collinear_3<K_base>                    Collinear_3;

  typedef Static_filters_predicates::Angle_3<K_base>                        Angle_3;
#endif // NOT CGAL_NO_ANGLE_3_STATIC_FILTERS
#ifndef CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS
  typedef Static_filters_predicates::Do_intersect_3<K_base,Self>            Do_intersect_3;
#endif // NOT CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS
  typedef Static_filters_predicates::Side_of_oriented_circle_2<K_base>      Side_of_oriented_circle_2;
  typedef Static_filters_predicates::Side_of_oriented_sphere_3<K_base>      Side_of_oriented_sphere_3;
  typedef Static_filters_predicates::Compare_squared_radius_3<K_base>       Compare_squared_radius_3;

  typedef Static_filters_predicates::Compare_weighted_squared_radius_3<K_base>     Compare_weighted_squared_radius_3;
  typedef Static_filters_predicates::Power_side_of_oriented_power_sphere_3<K_base>                          Power_side_of_oriented_power_sphere_3;

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(); }

  Orientation_3
  orientation_3_object() const
  { return Orientation_3(); }

  Collinear_3
  collinear_3_object() const
  { return Collinear_3(); }

#ifndef CGAL_NO_EQUAL_3_STATIC_FILTERS
 Equal_2
  equal_2_object() const
  { return Equal_2(); }

 Equal_3
  equal_3_object() const
  { return Equal_3(); }
#endif // NOT CGAL_NO_EQUAL_3_STATIC_FILTERS

#ifndef CGAL_NO_COMPARE_X_2_STATIC_FILTERS
 Compare_x_2
  compare_x_2_object() const
  { return Compare_x_2(); }

Compare_y_2
  compare_y_2_object() const
  { return Compare_y_2(); }
#endif // NOT CGAL_NO_COMPARE_Y_2_STATIC_FILTERS

#ifndef CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS
 Is_degenerate_3
  is_degenerate_3_object() const
  { return Is_degenerate_3(); }
#endif // NOT CGAL_NO_IS_DEGENERATE_3_STATIC_FILTERS

#ifndef CGAL_NO_ANGLE_3_STATIC_FILTERS
  Angle_3
  angle_3_object() const
  { return Angle_3(); }
#endif // NOT CGAL_NO_ANGLE_3_STATIC_FILTERS

#ifndef CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS
  Do_intersect_3
  do_intersect_3_object() const
  { return Do_intersect_3(); }
#endif // NOT CGAL_NO_DO_INTERSECT_3_STATIC_FILTERS

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }

  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const
  { return Side_of_oriented_sphere_3(); }

  Compare_squared_radius_3
  compare_squared_radius_3_object() const
  { return Compare_squared_radius_3(); }

  Power_side_of_oriented_power_sphere_3 power_side_of_oriented_power_sphere_3_object() const
  { return Power_side_of_oriented_power_sphere_3();}

  Compare_weighted_squared_radius_3
  compare_weighted_squared_radius_3_object() const
  { return Compare_weighted_squared_radius_3(); }


  enum { Has_static_filters = true };
}; // end of class template Static_filters<K_base, false>


// Here is the case when has_cheap_access_to_cartesian_coordinates is true,
// the default, used by Filtered_kernel<CK>.
// The K_base argument is supposed to provide exact primitives.
template < typename K_base>
class Static_filters<K_base, true> // has_cheap_access_to_cartesian_coordinates==true
  : public Static_filters<K_base, false>
{
  typedef Static_filters<K_base, true>              Self;

public:

  typedef Static_filters_predicates::Compare_y_at_x_2<K_base,Self>          Compare_y_at_x_2;

  // The following do not require filtering as they only do
  // comparisons.  We must be careful that *all* their function
  // operators must not do any comparisons.
  // In case we would like to avoid filtering only some of
  // the function operators, we need to make a new functors.
  typedef CommonKernelFunctors::Left_turn_2<Self>   Left_turn_2;

  typedef CartesianKernelFunctors::Less_x_2<Self>   Less_x_2;
  typedef CartesianKernelFunctors::Less_y_2<Self>   Less_y_2;
  typedef CartesianKernelFunctors::Less_xy_2<Self>  Less_xy_2;
  typedef CartesianKernelFunctors::Less_yx_2<Self>  Less_yx_2;

  typedef CartesianKernelFunctors::Less_x_3<Self>   Less_x_3;
  typedef CartesianKernelFunctors::Less_y_3<Self>   Less_y_3;
  typedef CartesianKernelFunctors::Less_z_3<Self>   Less_z_3;
  typedef CartesianKernelFunctors::Less_xy_3<Self>  Less_xy_3;
  typedef CartesianKernelFunctors::Less_xyz_3<Self> Less_xyz_3;

  typedef CartesianKernelFunctors::Compare_xy_2<Self>  Compare_xy_2;
  typedef CartesianKernelFunctors::Compare_x_3<Self>   Compare_x_3;
  typedef CartesianKernelFunctors::Compare_y_3<Self>   Compare_y_3;
  typedef CartesianKernelFunctors::Compare_z_3<Self>   Compare_z_3;
  typedef CartesianKernelFunctors::Compare_xy_3<Self>  Compare_xy_3;
  typedef CartesianKernelFunctors::Compare_xyz_3<Self> Compare_xyz_3;

  Compare_xy_2
  compare_xy_2_object() const
  { return Compare_xy_2(); }

  Compare_x_3
  compare_x_3_object() const
  { return Compare_x_3(); }

  Compare_y_3
  compare_y_3_object() const
  { return Compare_y_3(); }

  Compare_z_3
  compare_z_3_object() const
  { return Compare_z_3(); }

  Compare_xy_3
  compare_xy_3_object() const
  { return Compare_xy_3(); }

  Compare_xyz_3
  compare_xyz_3_object() const
  { return Compare_xyz_3(); }

  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(); }

  Less_x_2
  less_x_2_object() const
  { return Less_x_2(); }

  Less_y_2
  less_y_2_object() const
  { return Less_y_2(); }

  Less_xy_2
  less_xy_2_object() const
  { return Less_xy_2(); }

  Less_yx_2
  less_yx_2_object() const
  { return Less_yx_2(); }

  Less_x_3
  less_x_3_object() const
  { return Less_x_3(); }

  Less_y_3
  less_y_3_object() const
  { return Less_y_3(); }

  Less_z_3
  less_z_3_object() const
  { return Less_z_3(); }

  Less_xy_3
  less_xy_3_object() const
  { return Less_xy_3(); }

  Less_xyz_3
  less_xyz_3_object() const
  { return Less_xyz_3(); }

  Compare_y_at_x_2
  compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(); }

  // The two following are for degenerate cases, so I'll update them later.
  //
  // typedef Static_filters_predicates::Coplanar_orientation_3<Point_3, Orientation_2>
  //                                                   Coplanar_orientation_3;
  // typedef Static_filters_predicates::Side_of_bounded_circle_3<Point_3>
  //                                         Coplanar_side_of_bounded_circle_3;

  // Coplanar_orientation_3
  // coplanar_orientation_3_object() const
  // { return Coplanar_orientation_3(); }

  // Coplanar_side_of_bounded_circle_3
  // coplanar_side_of_bounded_circle_3_object() const
  // { return Coplanar_side_of_bounded_circle_3(); }
};

} } // namespace CGAL::internal

#endif // CGAL_INTERNAL_STATIC_FILTERS_H
