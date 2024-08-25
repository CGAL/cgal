// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel      <efif@post.tau.ac.il>
//            Eric Berberich <ericb@post.tau.ac.il>

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * A meradata traits-class decorator for the arrangement package. It counts the
 * number of invocations of traits-class functors. It is parameterized with
 * another traits class and inherits from it. For each traits method it
 * maintains a counter that counts the number of invocations into the method.
 *
 * It models all the concept that the original trais models.
 */

template <typename BaseTraits>
class Arr_counting_traits_2 : public BaseTraits {
public:
  enum Operation_id {
    COMPARE_X_2_OP = 0,
    COMPARE_XY_2_OP,
    CONSTRUCT_MIN_VERTEX_2_OP,
    CONSTRUCT_MAX_VERTEX_2_OP,
    IS_VERTICAL_2_OP,
    COMPARE_Y_AT_X_2_OP,
    EQUAL_2_POINTS_OP,
    EQUAL_2_CURVES_OP,
    COMPARE_Y_AT_X_LEFT_2_OP,
    COMPARE_Y_AT_X_RIGHT_2_OP,
    MAKE_X_MONOTONE_2_OP,
    SPLIT_2_OP,
    INTERSECT_2_OP,
    ARE_MERGEABLE_2_OP,
    MERGE_2_OP,
    CONSTRUCT_2_OPPOSITE_2_OP,
    COMPARE_ENDPOINTS_XY_2_OP,
    APPROXIMATE_2_COORD_OP,
    APPROXIMATE_2_POINT_OP,
    APPROXIMATE_2_CURVE_OP,
    PARAMETER_SPACE_IN_X_2_CURVE_END_OP,
    PARAMETER_SPACE_IN_X_2_POINT_OP,
    IS_ON_X_IDENTIFICATION_POINT_2_OP,
    IS_ON_X_IDENTIFICATION_CURVE_2_OP,
    COMPARE_Y_ON_BOUNDARY_2_OP,
    COMPARE_Y_NEAR_BOUNDARY_2_OP,
    PARAMETER_SPACE_IN_Y_2_CURVE_END_OP,
    PARAMETER_SPACE_IN_Y_2_POINT_OP,
    IS_ON_Y_IDENTIFICATION_2_POINT_OP,
    IS_ON_Y_IDENTIFICATION_2_CURVE_OP,
    COMPARE_X_ON_BOUNDARY_2_POINTS_OP,
    COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP,
    COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP,
    COMPARE_X_NEAR_BOUNDARY_2_OP,
    NUMBER_OF_OPERATIONS
  };

  /// \name Creation
  /// @{

  /*! Construct default */
  template <typename ... Args>
  Arr_counting_traits_2(Args ... args) : Base(std::forward<Args>(args)...) {}

  /*! Disable copy constructor.
   */
  Arr_counting_traits_2(const Arr_counting_traits_2&) = delete;

  /// @}

  /*! Obtain the counter of the given operation */
  std::size_t count(Operation_id id) const;

  /*! Print the compare_x counter */
  template <typename OutStream>
  OutStream& print(OutStream& os, Operation_id id) const;

  /// \name Types and functors inherited from the base
  /// @{

  using Has_left_category = typename Base::Has_left_category;
  using Has_merge_category = typename Base::Has_merge_category;
  using Has_do_intersect_category = typename Base::Has_do_intersect_category;

  using Left_side_category =
    typename internal::Arr_complete_left_side_category<Base>::Category;
  using Bottom_side_category =
    typename internal::Arr_complete_bottom_side_category<Base>::Category;
  using Top_side_category =
    typename internal::Arr_complete_top_side_category<Base>::Category;
  using Right_side_category =
    typename internal::Arr_complete_right_side_category<Base>::Category;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Curve_2 = typename Base::Curve_2;

  /// @}

  /// \name Obtain the appropriate functor
  /// @{

  Compare_x_2 compare_x_2_object() const;
  Compare_xy_2 compare_xy_2_object() const;
  Construct_min_vertex_2 construct_min_vertex_2_object() const;
  Construct_max_vertex_2 construct_max_vertex_2_object() const;
  Is_vertical_2 is_vertical_2_object() const;
  Compare_y_at_x_2 compare_y_at_x_2_object() const;
  Equal_2 equal_2_object() const;
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const;
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const;
  Make_x_monotone_2 make_x_monotone_2_object() const;
  Split_2 split_2_object() const;
  Intersect_2 intersect_2_object() const;
  Are_mergeable_2 are_mergeable_2_object() const;
  Merge_2 merge_2_object() const;
  Construct_opposite_2 construct_opposite_2_object() const;
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const;
  Approximate_2 approximate_2_object() const;
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const;
  Is_on_x_identification_2 is_on_x_identification_2_object() const;
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const;
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const;
  Is_on_y_identification_2 is_on_y_identification_2_object() const;
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const;

  /// @}

  /*! Clean all operation counters */
  void clear_counters();
};

template <typename OutStream, class BaseTraits>
inline OutStream& operator<<(OutStream& os,
                             const Arr_counting_traits_2<BaseTraits>& traits);
} //namespace CGAL
