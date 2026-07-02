/// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel    <efif@post.tau.ac.il>

namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * A metadata traits-class decorator for the arrangement package. It traces the
 * invocations of traits-class functors. It is parameterized with another traits
 * class and inherits from it. For each traits method it prints out its input
 * parameters and its output result
 *
 * It models all the concepts that the original traits models.
 */
template <typename BaseTraits>
class Arr_tracing_traits_2 {
public:
  enum Operation_id {
    COMPARE_X_2_OP = 0,
    COMPARE_XY_2_OP,
    CONSTRUCT_MIN_VERTEX_2_OP,
    CONSTRUCT_MAX_VERTEX_2_OP,
    IS_VERTICAL_2_OP,
    COMPARE_Y_AT_X_2_OP,
    EQUAL_POINTS_2_OP,
    EQUAL_CURVES_2_OP,
    COMPARE_Y_AT_X_LEFT_2_OP,
    COMPARE_Y_AT_X_RIGHT_2_OP,
    MAKE_X_MONOTONE_2_OP,
    SPLIT_2_OP,
    DO_INTERSECT_2_OP,
    INTERSECT_2_OP,
    ARE_MERGEABLE_2_OP,
    MERGE_2_OP,
    CONSTRUCT_2_OPPOSITE_2_OP,
    CONSTRUCT_POINT_2_OP,
    CONSTRUCT_POINT_2_XY_OP,
    CONSTRUCT_X_MONOTONE_CURVE_2_OP,
    CONSTRUCT_CURVE_2_OP,
    COMPARE_ENDPOINTS_XY_2_OP,
    APPROXIMATE_2_OP,
    PARAMETER_SPACE_IN_X_2_OP,
    IS_ON_X_IDENTIFICATION_2_OP,
    COMPARE_Y_ON_BOUNDARY_2_OP,
    COMPARE_Y_NEAR_BOUNDARY_2_OP,
    PARAMETER_SPACE_IN_Y_2_OP,
    IS_ON_Y_IDENTIFICATION_2_OP,
    COMPARE_X_ON_BOUNDARY_2_OP,
    COMPARE_X_NEAR_BOUNDARY_2_OP,
    NUMBER_OF_OPERATIONS
  };

public:
  /// \name Creation
  /// @{

  /*! constructs default.
   */
  template<typename ... Args>
  Arr_tracing_traits_2(Args ... args);

  /*! constructs from a shared pointer.
   * \param[in] traits the taits being traced.
   */
  Arr_tracing_traits_2(std::shared_ptr<BaseTraits> traits);

  /*! disables copy constructor.
   */
  Arr_tracing_traits_2(const Arr_tracing_traits_2&) = delete;

  /// @}

  /*! enables the trace of a traits operation
   * \param id the operation identifier
   */
  void enable_trace(Operation_id id);

  /*! enables the trace of all traits operations
   */
  void enable_all_traces();

  /*! disables the trace of a traits operation
   * \param id the operation identifier
   */
  void disable_trace(Operation_id id);

  /*! disables the trace of all traits operations
   */
  void disable_all_traces();

  /*! obtains a const reference to the traits being traced.
   */
  const Base& traits() const { return *m_base_traits; }

  /*! obtains a reference to the traits being traced.
   */
  Base& traits() { return *m_base_traits; }

  /*! obtains the smart pointer to the traits being traced.
   */
  std::shared_ptr<BaseTraits> shared_traits() const { return m_base_traits; }

  /// \name Types inherited from `BaseTraits`
  /// @{

  using Has_left_category = typename Base::Has_left_category;
  using Has_merge_category = typename Base::Has_merge_category;

  using Left_side_category = typename internal::Arr_complete_left_side_category<Base>::Category;
  using Bottom_side_category = typename internal::Arr_complete_bottom_side_category<Base>::Category;
  using Top_side_category = typename internal::Arr_complete_top_side_category<Base>::Category;
  using Right_side_category = typename internal::Arr_complete_right_side_category<Base>::Category;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  //! Defined only if the traits being traced models the concept `AosTraits_2`.
  using Curve_2 = typename Base::Curve_2;

  //! Defined only if the traits being traced models the concept `AosXMonotoneTraits_2`.
  using Multiplicity = typename Base::Multiplicity;

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
  Do_intersect_2 do_intersect_2_object() const;

  //! Supported only if the traits being traced models the concept `AosXMonotoneTraits_2`.
  Split_2 split_2_object() const;
  Intersect_2 intersect_2_object() const;
  Are_mergeable_2 are_mergeable_2_object() const;
  Merge_2 merge_2_object() const;

  //! Supported only if the traits being traced models the concept `AosTraits_2`.
  Make_x_monotone_2 make_x_monotone_2_object() const;

  //! Supported only if the traits being traced models the concept `AosDirectionalXMonotoneTraits_2`.
  Construct_opposite_2 construct_opposite_2_object() const;
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const;

  //! Supported only if the traits being traced models the concept `AosApproximateTraits_2`.
  Approximate_2 approximate_2_object() const;

  //! Supported only if the traits being traced models the concept `AosConstructPointTraits_2`
  Construct_point_2 construct_point_2_object() const;

  //! Supported only if the traits being traced models the concept `AosConstructXMonotoneCurveTraits_2`
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

  //! Supported only if the traits being traced models the concept `AosConstructCurveTraits_2`
  Construct_curve_2 construct_curve_2_object() const;

  /*! Supported only if the traits being traced models the concepts
   * `AosOpenBoundaryTraits_2` or `AosSphericalBoundaryTraits_2`.
   */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const;
  Is_on_x_identification_2 is_on_x_identification_2_object() const;
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const;
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const;
  Is_on_y_identification_2 is_on_y_identification_2_object() const;
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const;

  /// @}
};

template <typename OutputStream>
OutputStream& operator<<(OutputStream& os, Comparison_result cr);

} // namespace CGAL
