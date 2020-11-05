/*!
  \ingroup PkgKernel23Concepts
  \cgalConcept

  The concept of a <I>kernel</I> is defined by a set of requirements on
  the provision of certain types and access member functions to create
  objects of these types. The types are function object classes to be used
  within the algorithms and data structures of \cgal.
  This allows you to use any model of a kernel as a traits class in
  the \cgal algorithms and data structures, unless they require types
  beyond those provided by a kernel.

  A kernel provides types, construction objects, and generalized predicates.
  The former replace constructors of the kernel classes and constructive
  procedures in the kernel. There are also function objects replacing operators,
  especially for equality testing.

  \cgalHeading{Naming convention of constructions}
  All constructions which result type is a geometric object are prefixed by
  `Construct_`. If the result type is a number type, the name is prefixed by `Compute_`.
  When the result type is not determined, no prefix is used.

  \cgalHasModel `CGAL::Cartesian<FieldNumberType>`
  \cgalHasModel `CGAL::Homogeneous<RingNumberType>`
  \cgalHasModel `CGAL::Simple_cartesian<FieldNumberType>`
  \cgalHasModel `CGAL::Simple_homogeneous<RingNumberType>`
  \cgalHasModel `CGAL::Filtered_kernel<CK>`
  \cgalHasModel `CGAL::Exact_predicates_exact_constructions_kernel`
  \cgalHasModel `CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt`
  \cgalHasModel `CGAL::Exact_predicates_inexact_constructions_kernel`

  \sa `Kernel_d`
  \sa `CGAL::Ambient_dimension`
  \sa `CGAL::Feature_dimension`
  \sa `CGAL::Kernel_traits`

  \todo `Kernel::ConstructUnitNormal_3` as no model in the concept
  \todo `Kernel::CompareSquaredRadius_3` as no model in the concept
*/
class Kernel {
public:
  /// \name Types
  /// The following types describe the return types of predicates. They
  /// typically map to `bool` and \cgal kernel enum types, except when
  /// an interval arithmetic number type is used such as within the
  /// filtering kernels, in which case it is `Uncertain<bool>` or
  /// similar.
  /// @{

  /*!
    a model of `FieldNumberType`
  */
  typedef unspecified_type FT;

  /*!
    a model of `RingNumberType`
  */
  typedef unspecified_type RT;

  /*!
    `bool` or `Uncertain<bool>`
  */
  typedef unspecified_type Boolean;

  /*!
    `CGAL::Sign` or `Uncertain<CGAL::Sign>`
  */
  typedef unspecified_type Sign;

  /*!
    `CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>`
  */
  typedef unspecified_type Comparison_result;

  /*!
    `CGAL::Orientation` or `Uncertain<CGAL::Orientation>`
  */
  typedef unspecified_type Orientation;

  /*!
    `CGAL::Oriented_side` or `Uncertain<CGAL::Oriented_side>`
  */
  typedef unspecified_type Oriented_side;

  /*!
    `CGAL::Bounded_side` or `Uncertain<CGAL::Bounded_side>`
  */
  typedef unspecified_type Bounded_side;

  /*!
    `CGAL::Angle` or `Uncertain<CGAL::Angle>`
  */
  typedef unspecified_type Angle;

  /// @}

  /// \name Constants
  /// @{

  /*!
    A Boolean value indicating whether the predicates are filtered (as in
    `CGAL::Filtered_kernel`).  This helps propagating such decisions to traits
    classes which are built on top of a kernel, so that they can decide to
    filter their own predicates or not.
  */
  static const bool Has_filtered_predicates;

  /// @}

  /// \name Two-dimensional Coordinate Access
  /// @{

  /*!
    a model of `Kernel::CartesianConstIterator_2`
  */
  typedef unspecified_type Cartesian_const_iterator_2;

  /// @}

  /// \name Two-dimensional Geometric Objects
  /// @{

  /*!
    a model of `Kernel::Point_2`
  */
  typedef unspecified_type Point_2;

  /*!
    a model of `Kernel::WeightedPoint_2`
  */
  typedef unspecified_type Weighted_point_2;

  /*!
    a model of `Kernel::Vector_2`
  */
  typedef unspecified_type Vector_2;

  /*!
    a model of `Kernel::Direction_2`
  */
  typedef unspecified_type Direction_2;

  /*!
    a model of `Kernel::Line_2`
  */
  typedef unspecified_type Line_2;

  /*!
    a model of `Kernel::Ray_2`
  */
  typedef unspecified_type Ray_2;

  /*!
    a model of `Kernel::Segment_2`
  */
  typedef unspecified_type Segment_2;

  /*!
    a model of `Kernel::Triangle_2`
  */
  typedef unspecified_type Triangle_2;

  /*!
    a model of `Kernel::IsoRectangle_2`
  */
  typedef unspecified_type Iso_rectangle_2;

  /*!
    a model of `Kernel::Circle_2`
  */
  typedef unspecified_type Circle_2;

  /*!
    a model of `Kernel::Object_2`
  */
  typedef unspecified_type Object_2;

  /// @}

  /// \name Two-dimensional Constructions
  /// @{

  /*!
    a model of `Kernel::ConstructPoint_2`
  */
  typedef unspecified_type Construct_point_2;

  /*!
    a model of `Kernel::ConstructWeightedPoint_2`
  */
  typedef unspecified_type Construct_weighted_point_2;

  /*!
    a model of `Kernel::ConstructVector_2`
  */
  typedef unspecified_type Construct_vector_2;

  /*!
    a model of `Kernel::ConstructDirection_2`
  */
  typedef unspecified_type Construct_direction_2;

  /*!
    a model of `Kernel::ConstructSegment_2`
  */
  typedef unspecified_type Construct_segment_2;

  /*!
    a model of `Kernel::ConstructLine_2`
  */
  typedef unspecified_type Construct_line_2;

  /*!
    a model of `Kernel::ConstructRay_2`
  */
  typedef unspecified_type Construct_ray_2;

  /*!
    a model of `Kernel::ConstructCircle_2`
  */
  typedef unspecified_type Construct_circle_2;

  /*!
    a model of `Kernel::ConstructTriangle_2`
  */
  typedef unspecified_type Construct_triangle_2;

  /*!
    a model of `Kernel::ConstructIsoRectangle_2`
  */
  typedef unspecified_type Construct_iso_rectangle_2;

  /*!
    a model of `Kernel::ConstructObject_2`
  */
  typedef unspecified_type Construct_object_2;

  /*!
    a model of `Kernel::ConstructScaledVector_2`
  */
  typedef unspecified_type Construct_scaled_vector_2;

  /*!
    a model of `Kernel::ConstructDifferenceOfVectors_2`
  */
  typedef unspecified_type Construct_difference_of_vectors_2;

  /*!
    a model of `Kernel::ConstructSumOfVectors_2`
  */
  typedef unspecified_type Construct_sum_of_vectors_2;

  /*!
    a model of `Kernel::ConstructDividedVector_2`
  */
  typedef unspecified_type Construct_divided_vector_2;

  /*!
    a model of `Kernel::ConstructTranslatedPoint_2`
  */
  typedef unspecified_type Construct_translated_point_2;

  /*!
    a model of `Kernel::ConstructPointOn_2`
  */
  typedef unspecified_type Construct_point_on_2;

  /*!
    a model of `Kernel::ConstructProjectedPoint_2`
  */
  typedef unspecified_type Construct_projected_point_2;

  /*!
    a model of `Kernel::ConstructProjectedXYPoint_2`
  */
  typedef unspecified_type Construct_projected_xy_point_2;

  /*!
    a model of `Kernel::ConstructSecondPoint_2`
    \todo This is used in CGAL::Ray_2<K>::second_point() but is not documented
  */
  typedef unspecified_type Construct_second_point_2;

  /*!
    a model of `Kernel::ConstructCartesianConstIterator_2`
  */
  typedef unspecified_type Construct_cartesian_const_iterator_2;

  /*!
    a model of `Kernel::ConstructVertex_2`
  */
  typedef unspecified_type Construct_vertex_2;

  /*!
    a model of `Kernel::ConstructBbox_2`
  */
  typedef unspecified_type Construct_bbox_2;

  /*!
    a model of `Kernel::ConstructPerpendicularVector_2`
  */
  typedef unspecified_type Construct_perpendicular_vector_2;

  /*!
    a model of `Kernel::ConstructPerpendicularDirection_2`
  */
  typedef unspecified_type Construct_perpendicular_direction_2;

  /*!
    a model of `Kernel::ConstructPerpendicularLine_2`
  */
  typedef unspecified_type Construct_perpendicular_line_2;


  /*!
    a model of `Kernel::ConstructMidpoint_2`
  */
  typedef unspecified_type Construct_midpoint_2;

  /*!
    a model of `Kernel::ConstructEquidistantLine_3`
  */
  typedef unspecified_type Construct_equidistant_line_3;

  /*!
    a model of `Kernel::ConstructCenter_2`
  */
  typedef unspecified_type Construct_center_2;

  /*!
    a model of `Kernel::ConstructCentroid_2`
  */
  typedef unspecified_type Construct_centroid_2;

  /*!
    a model of `Kernel::ConstructBarycenter_2`
  */
  typedef unspecified_type Construct_barycenter_2;

  /*!
    a model of `Kernel::ConstructCircumcenter_2`
  */
  typedef unspecified_type Construct_circumcenter_2;

  /*!
    a model of `Kernel::ConstructWeightedCircumcenter_2`
  */
  typedef unspecified_type Construct_weighted_circumcenter_2;

  /*!
    a model of `Kernel::ConstructBisector_2`
  */
  typedef unspecified_type Construct_bisector_2;

  /*!
    a model of `Kernel::ConstructRadicalAxis_2`
  */
  typedef unspecified_type Construct_radical_axis_2;

  /*!
    a model of `Kernel::ConstructRadicalLine_2`
  */
  typedef unspecified_type Construct_radical_line_2;

  /*!
    a model of `Kernel::ConstructOppositeDirection_2`
  */
  typedef unspecified_type Construct_opposite_direction_2;

  /*!
    a model of `Kernel::ConstructOppositeSegment_2`
  */
  typedef unspecified_type Construct_opposite_segment_2;

  /*!
    a model of `Kernel::ConstructOppositeRay_2`
  */
  typedef unspecified_type Construct_opposite_ray_2;

  /*!
    a model of `Kernel::ConstructOppositeLine_2`
  */
  typedef unspecified_type Construct_opposite_line_2;

  /*!
    a model of `Kernel::ConstructOppositeTriangle_2`
  */
  typedef unspecified_type Construct_opposite_triangle_2;

  /*!
    a model of `Kernel::ConstructOppositeCircle_2`
  */
  typedef unspecified_type Construct_opposite_circle_2;

  /*!
    a model of `Kernel::ConstructOppositeVector_2`
  */
  typedef unspecified_type Construct_opposite_vector_2;

  /*!
    a model of `Kernel::ConstructSource_2`
  */
  typedef unspecified_type  Construct_source_2;

  /*!
    a model of `Kernel::ConstructTarget_2`
  */
  typedef unspecified_type  Construct_target_2;

  /*!
    a model of `Kernel::ConstructMaxVertex_2`
  */
  typedef unspecified_type  Construct_max_vertex_2;

  /*!
    a model of `Kernel::ConstructMinVertex_2`
  */
  typedef unspecified_type  Construct_min_vertex_2;

  /*!
    a model of `Kernel::Intersect_2`
  */
  typedef unspecified_type Intersect_2;

  /*!
  a model of `Kernel::Assign_2`
  */
  typedef unspecified_type Assign_2;

  /*!
    a model of `Kernel::ComputeWeight_2`
  */
  typedef unspecified_type Compute_weight_2;

  /*!
    a model of `Kernel::ComputeX_2`
  */
  typedef unspecified_type Compute_x_2;

  /*!
    a model of `Kernel::ComputeY_2`
  */
  typedef unspecified_type Compute_y_2;

  /*!
    a model of `Kernel::ComputeA_2`
  */
  typedef unspecified_type  Compute_a_2;

  /*!
    a model of `Kernel::ComputeB_2`
  */
  typedef unspecified_type  Compute_b_2;

  /*!
    a model of `Kernel::ComputeC_2`
  */
  typedef unspecified_type  ComputeC_2;

  /*!
    a model of `Kernel::ComputeDx_2`
  */
  typedef unspecified_type  Compute_dx_2;

  /*!
    a model of `Kernel::ComputeDy_2`
  */
  typedef unspecified_type  Compute_dy_2;

  /*!
    a model of `Kernel::ComputeHx_2`
  */
  typedef unspecified_type  Compute_hx_2;

  /*!
    a model of `Kernel::ComputeHy_2`
  */
  typedef unspecified_type  Compute_hy_2;

  /*!
    a model of `Kernel::ComputeHw_2`
  */
  typedef unspecified_type  Compute_hw_2;

  /*!
    a model of `Kernel::ComputeLInfinityDistance_2`
  */
  typedef unspecified_type  Compute_L_infinity_distance_2;

  /*!
    a model of `Kernel::ComputeXmax_2`
  */
  typedef unspecified_type  Compute_xmax_2;

  /*!
    a model of `Kernel::ComputeXmin_2`
  */
  typedef unspecified_type  Compute_xmin_2;

  /*!
    a model of `Kernel::ComputeYmax_2`
  */
  typedef unspecified_type  Compute_ymax_2;

  /*!
    a model of `Kernel::ComputeYmin_2`
  */
  typedef unspecified_type  Compute_ymin_2;

  /*!
    a model of `Kernel::ComputeYAtX_2`
  */
  typedef unspecified_type  Compute_y_at_x_2;

  /*!
    a model of `Kernel::ComputeSquaredDistance_2`
  */
  typedef unspecified_type Compute_squared_distance_2;

  /*!
    a model of `Kernel::ComputeSquaredLength_2`
  */
  typedef unspecified_type Compute_squared_length_2;

  /*!
    a model of `Kernel::ComputeSquaredRadius_2`
  */
  typedef unspecified_type Compute_squared_radius_2;

  /*!
    a model of `Kernel::ComputeSquaredRadiusSmallestOrthogonalCircle_2`
  */
  typedef unspecified_type Compute_squared_radius_smallest_orthogonal_circle_2;

  /*!
    a model of `Kernel::ComputePowerProduct_2`
  */
  typedef unspecified_type Compute_power_product_2;

  /*!
    a model of `Kernel::ComputeArea_2`
  */
  typedef unspecified_type Compute_area_2;

  /*!
    a model of `Kernel::ComputeDeterminant_2`
  */
  typedef unspecified_type Compute_determinant_2;

  /*!
    a model of `Kernel::ComputeScalarProduct_2`
  */
  typedef unspecified_type Compute_scalar_product_2;

  /// @}

  /// \name Two-dimensional Generalized Predicates
  /// @{

  /*!
    a model of `Kernel::Angle_2`
  */
  typedef unspecified_type Angle_2;

  /*!
    a model of `Kernel::Equal_2`
  */
  typedef unspecified_type Equal_2;

  /*!
    a model of `Kernel::EqualX_2`
  */
  typedef unspecified_type Equal_x_2;

  /*!
    a model of `Kernel::EqualY_2`
  */
  typedef unspecified_type Equal_y_2;

  /*!
    a model of `Kernel::LessX_2`
  */
  typedef unspecified_type Less_x_2;

  /*!
    a model of `Kernel::LessY_2`
  */
  typedef unspecified_type Less_y_2;

  /*!
    a model of `Kernel::LessXY_2`
  */
  typedef unspecified_type Less_xy_2;

  /*!
    a model of `Kernel::LessYX_2`
  */
  typedef unspecified_type Less_yx_2;

  /*!
    a model of `Kernel::CompareX_2`
  */
  typedef unspecified_type Compare_x_2;

  /*!
    a model of `Kernel::CompareXAtY_2`
  */
  typedef unspecified_type Compare_x_at_y_2;

  /*!
    a model of `Kernel::CompareY_2`
  */
  typedef unspecified_type Compare_y_2;

  /*!
    a model of `Kernel::CompareXY_2`
  */
  typedef unspecified_type Compare_xy_2;

  /*!
    a model of `Kernel::CompareYX_2`
  */
  typedef unspecified_type Compare_yx_2;

  /*!
    a model of `Kernel::CompareYAtX_2`
  */
  typedef unspecified_type Compare_y_at_x_2;

  /*!
    a model of `Kernel::CompareDistance_2`
  */
  typedef unspecified_type Compare_distance_2;

  /*!
    a model of `Kernel::ComparePowerDistance_2`
  */
  typedef unspecified_type Compare_power_distance_2;

  /*!
    a model of `Kernel::CompareSquaredDistance_2`
  */
  typedef unspecified_type Compare_square_distance_2;

  /*!
    a model of `Kernel::CompareAngleWithXAxis_2`
  */
  typedef unspecified_type Compare_angle_with_x_axis_2;

  /*!
    a model of `Kernel::CompareSignedDistanceToLine_2`
  */
  typedef unspecified_type Compare_signed_distance_to_line_2;

  /*!
    a model of `Kernel::CompareSlope_2`
  */
  typedef unspecified_type Compare_slope_2;

  /*!
    a model of `Kernel::LessDistanceToPoint_2`
  */
  typedef unspecified_type Less_distance_to_point_2;

  /*!
    a model of `Kernel::LessSignedDistanceToLine_2`
  */
  typedef unspecified_type Less_signed_distance_to_line_2;

  /*!
    a model of `Kernel::LessRotateCCW_2`
  */
  typedef unspecified_type Less_rotate_ccw_2;

  /*!
    a model of `Kernel::LeftTurn_2`
  */
  typedef unspecified_type Left_turn_2;

  /*!
    a model of `Kernel::Collinear_2`
  */
  typedef unspecified_type Collinear_2;

  /*!
    a model of `Kernel::Orientation_2`
  */
  typedef unspecified_type Orientation_2;

  /*!
    a model of `Kernel::SideOfOrientedCircle_2`
  */
  typedef unspecified_type Side_of_oriented_circle_2;

  /*!
    a model of `Kernel::SideOfBoundedCircle_2`
  */
  typedef unspecified_type Side_of_bounded_circle_2;

  /*!
    a model of `Kernel::PowerSideOfOrientedPowerCircle_2`
  */
  typedef unspecified_type Power_side_of_oriented_power_circle_2;

  /*!
    a model of `Kernel::PowerSideOfBoundedPowerCircle_2`
  */
  typedef unspecified_type Power_side_of_bounded_power_circle_2;

  /*!
    a model of `Kernel::IsHorizontal_2`
  */
  typedef unspecified_type Is_horizontal_2;

  /*!
    a model of `Kernel::IsVertical_2`
  */
  typedef unspecified_type Is_vertical_2;

  /*!
    a model of `Kernel::IsDegenerate_2`
  */
  typedef unspecified_type Is_degenerate_2;

  /*!
    a model of `Kernel::HasOn_2`
  */
  typedef unspecified_type Has_on_2;

  /*!
    a model of `Kernel::CollinearHasOn_2`
  */
  typedef unspecified_type Collinear_has_on_2;

  /*!
    a model of `Kernel::HasOnBoundedSide_2`
  */
  typedef unspecified_type Has_on_bounded_side_2;

  /*!
    a model of `Kernel::HasOnUnboundedSide_2`
  */
  typedef unspecified_type Has_on_unbounded_side_2;

  /*!
    a model of `Kernel::HasOnBoundary_2`
  */
  typedef unspecified_type Has_on_boundary_2;

  /*!
    a model of `Kernel::HasOnPositiveSide_2`
  */
  typedef unspecified_type Has_on_positive_side_2;

  /*!
    a model of `Kernel::HasOnNegativeSide_2`
  */
  typedef unspecified_type Has_on_negative_side_2;

  /*!
    a model of `Kernel::OrientedSide_2`
  */
  typedef unspecified_type Oriented_side_2;

  /*!
    a model of `Kernel::BoundedSide_2`
  */
  typedef unspecified_type Bounded_side_2;

  /*!
    a model of `Kernel::AreParallel_2`
  */
  typedef unspecified_type Are_parallel_2 ;

  /*!
    a model of `Kernel::AreOrderedAlongLine_2`
  */
  typedef unspecified_type Are_ordered_along_line_2 ;

  /*!
    a model of `Kernel::AreStrictlyOrderedAlongLine_2`
  */
  typedef unspecified_type Are_strictly_ordered_along_line_2;

  /*!
    a model of `Kernel::CollinearAreOrderedAlongLine_2`
  */
  typedef unspecified_type Collinear_are_ordered_along_line_2;

  /*!
    a model of `Kernel::CollinearAreStrictlyOrderedAlongLine_2`
  */
  typedef unspecified_type Collinear_are_strictly_ordered_along_line_2;

  /*!
    a model of `Kernel::CounterclockwiseInBetween_2`
  */
  typedef unspecified_type Counterclockwise_in_between_2;

  /*!
    a model of `Kernel::DoIntersect_2`
  */
  typedef unspecified_type Do_intersect_2;

  /// @}

  /// \name Three-dimensional Coordinate Access
  /// @{

  /*!
    a model of `Kernel::CartesianConstIterator_3`
  */
  typedef unspecified_type Cartesian_const_iterator_3;

  /// @}

  /// \name Three-dimensional Geometric Objects
  /// @{

  /*!
    a model of `Kernel::Point_3`
  */
  typedef unspecified_type Point_3;

  /*!
    a model of `Kernel::WeightedPoint_3`
  */
  typedef unspecified_type Weighted_point_3;

  /*!
    a model of `Kernel::Vector_3`
  */
  typedef unspecified_type Vector_3;

  /*!
    a model of `Kernel::Direction_3`
  */
  typedef unspecified_type Direction_3;

  /*!
    a model of `Kernel::IsoCuboid_3`
  */
  typedef unspecified_type Iso_cuboid_3;

  /*!
    a model of `Kernel::Line_3`
  */
  typedef unspecified_type Line_3;

  /*!
    a model of `Kernel::Ray_3`
  */
  typedef unspecified_type Ray_3;

  /*!
    a model of `Kernel::Circle_3`
  */
  typedef unspecified_type Circle_3;

  /*!
    a model of `Kernel::Sphere_3`
  */
  typedef unspecified_type Sphere_3;

  /*!
    a model of `Kernel::Segment_3`
  */
  typedef unspecified_type Segment_3;

  /*!
    a model of `Kernel::Plane_3`
  */
  typedef unspecified_type Plane_3;

  /*!
    a model of `Kernel::Triangle_3`
  */
  typedef unspecified_type Triangle_3;

  /*!
    a model of `Kernel::Tetrahedron_3`
  */
  typedef unspecified_type Tetrahedron_3;

  /*!
    a model of `Kernel::Object_3`
  */
  typedef unspecified_type Object_3;

  /// @}

  /// \name Three-dimensional Constructions
  /// @{

  /*!
    a model of `Kernel::ConstructPoint_3`
  */
  typedef unspecified_type Construct_point_3;

  /*!
    a model of `Kernel::ConstructWeightedPoint_3`
  */
  typedef unspecified_type Construct_weighted_point_3;

  /*!
    a model of `Kernel::ConstructVector_3`
  */
  typedef unspecified_type Construct_vector_3;

  /*!
    a model of `Kernel::ConstructDirection_3`
  */
  typedef unspecified_type Construct_direction_3;

  /*!
    a model of `Kernel::ConstructPlane_3`
  */
  typedef unspecified_type Construct_plane_3;

  /*!
    a model of `Kernel::ConstructIsoCuboid_3`
  */
  typedef unspecified_type Construct_iso_cuboid_3;

  /*!
    a model of `Kernel::ConstructLine_3`
  */
  typedef unspecified_type Construct_line_3;

  /*!
    a model of `Kernel::ConstructRay_3`
  */
  typedef unspecified_type Construct_ray_3;

  /*!
    a model of `Kernel::ConstructSphere_3`
  */
  typedef unspecified_type Construct_sphere_3;

  /*!
    a model of `Kernel::ConstructSegment_3`
  */
  typedef unspecified_type Construct_segment_3;

  /*!
    a model of `Kernel::ConstructTriangle_3`
  */
  typedef unspecified_type Construct_triangle_3;

  /*!
    a model of `Kernel::ConstructTetrahedron_3`
  */
  typedef unspecified_type Construct_tetrahedron_3;

  /*!
    a model of `Kernel::ConstructCircle_3`
  */
  typedef unspecified_type Construct_circle_3;

  /*!
    a model of `Kernel::ConstructObject_3`
  */
  typedef unspecified_type Construct_object_3;

  /*!
    a model of `Kernel::ConstructScaledVector_3`
  */
  typedef unspecified_type Construct_scaled_vector_3;

  /*!
    a model of `Kernel::ConstructTranslatedPoint_3`
  */
  typedef unspecified_type Construct_translated_point_3;

  /*!
    a model of `Kernel::ConstructPointOn_3`
  */
  typedef unspecified_type Construct_point_on_3;

  /*!
    a model of `Kernel::ConstructProjectedPoint_3`
  */
  typedef unspecified_type Construct_projected_point_3;

  /*!
    a model of `Kernel::ConstructLiftedPoint_3`
  */
  typedef unspecified_type Construct_lifted_point_3;

  /*!
    a model of `Kernel::ConstructCartesianConstIterator_3`
  */
  typedef unspecified_type Construct_cartesian_const_iterator_3;

  /*!
    a model of `Kernel::ConstructVertex_3`
  */
  typedef unspecified_type Construct_vertex_3;

  /*!
    a model of `Kernel::ConstructSource_3`
  */
  typedef unspecified_type  Construct_source_3;

  /*!
    a model of `Kernel::ConstructTarget_3`
  */
  typedef unspecified_type  Construct_target_3;

  /*!
    a model of `Kernel::ConstructMaxVertex_3`
  */
  typedef unspecified_type  Construct_max_vertex_3;

  /*!
    a model of `Kernel::ConstructMinVertex_3`
  */
  typedef unspecified_type  Construct_min_vertex_3;

  /*!
    a model of `Kernel::ConstructBbox_3`
  */
  typedef unspecified_type Construct_bbox_3;

  /*!
    a model of `Kernel::ConstructSupportingPlane_3`
  */
  typedef unspecified_type Construct_supporting_plane_3;

  /*!
    a model of `Kernel::ConstructOrthogonalVector_3`
  */
  typedef unspecified_type Construct_orthogonal_vector_3;

  /*!
    a model of `Kernel::ConstructBaseVector_3`
  */
  typedef unspecified_type Construct_base_vector_3;

  /*!
    a model of `Kernel::ConstructDifferenceOfVectors_3`
  */
  typedef unspecified_type Construct_difference_of_vectors_3;

  /*!
    a model of `Kernel::ConstructSumOfVectors_3`
  */
  typedef unspecified_type Construct_sum_of_vectors_3;

  /*!
    a model of `Kernel::ConstructDividedVector_3`
  */
  typedef unspecified_type Construct_divided_vector_3;

  /*!
    a model of `Kernel::ConstructNormal_3`
  */
  typedef unspecified_type Construct_normal_3;

  /*!
    a model of `Kernel::ConstructPerpendicularPlane_3`
  */
  typedef unspecified_type Construct_perpendicular_plane_3;

  /*!
    a model of `Kernel::ConstructRadicalPlane_3`
  */
  typedef unspecified_type Construct_radical_plane_3;

  /*!
    a model of `Kernel::ConstructPerpendicularLine_3`
  */
  typedef unspecified_type Construct_perpendicular_line_3;

  /*!
    a model of `Kernel::ConstructMidpoint_3`
  */
  typedef unspecified_type Construct_midpoint_3;

  /*!
    a model of `Kernel::ConstructCenter_3`
  */
  typedef unspecified_type Construct_center_3;

  /*!
    a model of `Kernel::ConstructCentroid_3`
  */
  typedef unspecified_type Construct_centroid_3;

  /*!
    a model of `Kernel::ConstructCircumcenter_3`
  */
  typedef unspecified_type Construct_circumcenter_3;

  /*!
    a model of `Kernel::ConstructWeightedCircumcenter_3`
  */
  typedef unspecified_type Construct_weighted_circumcenter_3;

  /*!
    a model of `Kernel::ConstructBisector_3`
  */
  typedef unspecified_type Construct_bisector_3;

  /*!
    a model of `Kernel::ConstructCrossProductVector_3`
  */
  typedef unspecified_type Construct_cross_product_vector_3;

  /*!
    a model of `Kernel::ConstructBarycenter_3`
  */
  typedef unspecified_type Construct_barycenter_3;

  /*!
    a model of `Kernel::ConstructSecondPoint_3`
    \todo This is used in CGAL::Ray_2<K>::second_point() but is not documented
  */
  typedef unspecified_type Construct_second_point_3;

  /*!
    a model of `Kernel::ConstructOppositeDirection_3`
  */
  typedef unspecified_type Construct_opposite_direction_3;

  /*!
    a model of `Kernel::ConstructOppositeSegment_3`
  */
  typedef unspecified_type Construct_opposite_segment_3;

  /*!
    a model of `Kernel::ConstructOppositeRay_3`
  */
  typedef unspecified_type Construct_opposite_ray_3;

  /*!
    a model of `Kernel::ConstructOppositeLine_3`
  */
  typedef unspecified_type Construct_opposite_line_3;

  /*!
    a model of `Kernel::ConstructOppositePlane_3`
  */
  typedef unspecified_type Construct_opposite_plane_3;

  /*!
    a model of `Kernel::ConstructOppositeSphere_3`
  */
  typedef unspecified_type Construct_opposite_sphere_3;

  /*!
    a model of `Kernel::ConstructOppositeVector_3`
  */
  typedef unspecified_type Construct_opposite_vector_3;

  /*!
    a model of `Kernel::Intersect_3`
  */
  typedef unspecified_type Intersect_3;

  /*!
    a model of `Kernel::Assign_3`
  */
  typedef unspecified_type Assign_3;

  /*!
    a model of `Kernel::ComputeWeight_3`
  */
  typedef unspecified_type Compute_weight_3;

  /*!
    a model of `Kernel::ComputeX_3`
  */
  typedef unspecified_type Compute_x_3;

  /*!
    a model of `Kernel::ComputeY_3`
  */
  typedef unspecified_type Compute_y_3;

  /*!
    a model of `Kernel::ComputeZ_3`
  */
  typedef unspecified_type Compute_z_3;

  /*!
    a model of `Kernel::ComputeA_3`
  */
  typedef unspecified_type  Compute_a_3;

  /*!
    a model of `Kernel::ComputeB_3`
  */
  typedef unspecified_type  Compute_b_3;

  /*!
    a model of `Kernel::ComputeC_3`
  */
  typedef unspecified_type  Compute_c_3;

  /*!
    a model of `Kernel::ComputeD_3`
  */
  typedef unspecified_type  Compute_d_3;

  /*!
    a model of `Kernel::ComputeHx_3`
  */
  typedef unspecified_type  Compute_hx_3;

  /*!
    a model of `Kernel::ComputeHy_3`
  */
  typedef unspecified_type  Compute_hy_3;

  /*!
    a model of `Kernel::ComputeHz_3`
  */
  typedef unspecified_type  Compute_hz_3;

  /*!
    a model of `Kernel::ComputeHw_3`
  */
  typedef unspecified_type  Compute_hw_3;

  /*!
    a model of `Kernel::ComputeLInfinityDistance_3`
  */
  typedef unspecified_type  Compute_L_infinity_distance_3;

  /*!
    a model of `Kernel::ComputeDx_3`
  */
  typedef unspecified_type  Compute_dx_3;

  /*!
    a model of `Kernel::ComputeDy_3`
  */
  typedef unspecified_type  Compute_dy_3  ;

  /*!
    a model of `Kernel::ComputeDz_3`
  */
  typedef unspecified_type  Compute_dz_3;

  /*!
    a model of `Kernel::ComputeXmax_3`
  */
  typedef unspecified_type  Compute_xmax_3;

  /*!
    a model of `Kernel::ComputeXmin_3`
  */
  typedef unspecified_type  Compute_xmin_3;

  /*!
    a model of `Kernel::ComputeYmax_3`
  */
  typedef unspecified_type  Compute_ymax_3;

  /*!
    a model of `Kernel::ComputeYmin_3`
  */
  typedef unspecified_type  Compute_ymin_3;

  /*!
    a model of `Kernel::ComputeZmax_3`
  */
  typedef unspecified_type  Compute_zmax_3;

  /*!
    a model of `Kernel::ComputeZmin_3`
  */
  typedef unspecified_type  Compute_zmin_3;

  /*!
    a model of `Kernel::ComputeArea_3`
  */
  typedef unspecified_type Compute_area_3;

  /*!
    a model of `Kernel::ComputeSquaredArea_3`
  */
  typedef unspecified_type Compute_squared_area_3;

  /*!
    a model of `Kernel::ComputeAreaDividedByPi_3`
  */
  typedef unspecified_type Compute_area_divided_by_pi_3;

  /*!
    a model of `Kernel::ComputeApproximateArea_3`
  */
  typedef unspecified_type Compute_approximate_area_3;

  /*!
    a model of `Kernel::ComputeApproximateAngle_3`
  */
  typedef unspecified_type Compute_approximate_angle_3;

  /*!
    a model of `Kernel::ComputeApproximateDihedralAngle_3`
  */
  typedef unspecified_type Compute_approximate_dihedral_angle_3;

  /*!
    a model of `Kernel::ComputeDeterminant_3`
  */
  typedef unspecified_type Compute_determinant_3;

  /*!
    a model of `Kernel::ComputeScalarProduct_3`
  */
  typedef unspecified_type Compute_scalar_product_3;

  /*!
    a model of `Kernel::ComputeSquaredDistance_3`
  */
  typedef unspecified_type Compute_squared_distance_3;

  /*!
    a model of `Kernel::ComputeSquaredLength_3`
  */
  typedef unspecified_type Compute_squared_length_3;

  /*!
    a model of `Kernel::ComputeSquaredLengthDividedByPiSquare_3`
  */
  typedef unspecified_type Compute_squared_length_divided_by_pi_square_3;

  /*!
    a model of `Kernel::ComputeApproximateSquaredLength_3`
  */
  typedef unspecified_type Compute_approximate_squared_length_3;

  /*!
    a model of `Kernel::ComputeSquaredRadius_3`
  */
  typedef unspecified_type Compute_squared_radius_3;

  /*!
    a model of `Kernel::ComputeVolume_3`
  */
  typedef unspecified_type Compute_volume_3;

  /*!
    a model of `Kernel::ComputePowerProduct_3`
  */
  typedef unspecified_type Compute_power_product_3;

  /*!
    a model of `Kernel::ComputePowerDistanceToPowerSphere_3`
  */
  typedef unspecified_type Compute_power_distance_to_power_sphere_3;

  /*!
    a model of `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3`
  */
  typedef unspecified_type Compute_squared_radius_smallest_orthogonal_sphere_3;

  /// @}

  /// \name Three-dimensional Generalized Predicates
  /// @{

  /*!
    a model of `Kernel::Angle_3`
  */
  typedef unspecified_type Angle_3;

  /*!
    a model of `Kernel::Equal_3`
  */
  typedef unspecified_type Equal_3;

  /*!
    a model of `Kernel::EqualX_3`
  */
  typedef unspecified_type Equal_x_3;

  /*!
    a model of `Kernel::EqualY_3`
  */
  typedef unspecified_type Equal_y_3;

  /*!
    a model of `Kernel::EqualZ_3`
  */
  typedef unspecified_type Equal_z_3;

  /*!
    a model of `Kernel::EqualXY_3`
  */
  typedef unspecified_type Equal_xy_3;

  /*!
    a model of `Kernel::LessX_3`
  */
  typedef unspecified_type Less_x_3;

  /*!
    a model of `Kernel::LessY_3`
  */
  typedef unspecified_type Less_y_3;

  /*!
    a model of `Kernel::LessZ_3`
  */
  typedef unspecified_type Less_z_3;

  /*!
    a model of `Kernel::LessXY_3`
  */
  typedef unspecified_type Less_xy_3;

  /*!
    a model of `Kernel::LessXYZ_3`
  */
  typedef unspecified_type Less_xyz_3;

  /*!
    a model of `Kernel::CompareWeightedSquaredRadius_3`
  */
  typedef unspecified_type Compare_weighted_squared_radius_3;

  /*!
    a model of `Kernel::CompareX_3`
  */
  typedef unspecified_type Compare_x_3;

  /*!
    a model of `Kernel::CompareY_3`
  */
  typedef unspecified_type Compare_y_3;

  /*!
    a model of `Kernel::CompareZ_3`
  */
  typedef unspecified_type Compare_z_3;

  /*!
    a model of `Kernel::CompareXY_3`
  */
  typedef unspecified_type Compare_xy_3;

  /*!
    a model of `Kernel::CompareXYZ_3`
  */
  typedef unspecified_type Compare_xyz_3;

  /*!
    a model of `Kernel::CompareSlope_3`
  */
  typedef unspecified_type Compare_slope_3;

  /*!
    a model of `Kernel::CompareSquaredDistance_3`
  */
  typedef unspecified_type Compare_squared_distance_3;

  /*!
    a model of `Kernel::LessSignedDistanceToPlane_3`
  */
  typedef unspecified_type Less_signed_distance_to_plane_3;

  /*!
    a model of `Kernel::LessDistanceToPoint_3`
  */
  typedef unspecified_type Less_distance_to_point_3;

  /*!
    a model of `Kernel::CompareDihedralAngle_3`
  */
  typedef unspecified_type Compare_dihedral_angle_3;

  /*!
    a model of `Kernel::CompareDistance_3`
  */
  typedef unspecified_type Compare_distance_3;

  /*!
    a model of `Kernel::ComparePowerDistance_3`
  */
  typedef unspecified_type Compare_power_distance_3;

  /*!
    a model of `Kernel::Collinear_3`
  */
  typedef unspecified_type Collinear_3;

  /*!
    a model of `Kernel::Coplanar_3`
  */
  typedef unspecified_type Coplanar_3;

  /*!
    a model of `Kernel::Orientation_3`
  */
  typedef unspecified_type Orientation_3;

  /*!
    a model of `Kernel::CoplanarOrientation_3`
  */
  typedef unspecified_type Coplanar_orientation_3;

  /*!
    a model of `Kernel::CoplanarSideOfBoundedCircle_3`
  */
  typedef unspecified_type Coplanar_side_of_bounded_circle_3;

  /*!
    a model of `Kernel::SideOfOrientedSphere_3`
  */
  typedef unspecified_type Side_of_oriented_sphere_3;

  /*!
    a model of `Kernel::SideOfBoundedSphere_3`
  */
  typedef unspecified_type Side_of_bounded_sphere_3;

  /*!
    a model of `Kernel::PowerSideOfOrientedPowerSphere_3`
  */
  typedef unspecified_type Power_side_of_oriented_power_sphere_3;

  /*!
    a model of `Kernel::PowerSideOfBoundedPowerSphere_3`
  */
  typedef unspecified_type Power_side_of_bounded_power_sphere_3;

  /*!
    a model of `Kernel::IsDegenerate_3`
  */
  typedef unspecified_type Is_degenerate_3;

  /*!
    a model of `Kernel::HasOn_3`
  */
  typedef unspecified_type Has_on_3;

  /*!
    a model of `Kernel::HasOnBoundedSide_3`
  */
  typedef unspecified_type Has_on_bounded_side_3;

  /*!
    a model of `Kernel::HasOnUnboundedSide_3`
  */
  typedef unspecified_type Has_on_unbounded_side_3;

  /*!
    a model of `Kernel::HasOnBoundary_3`
  */
  typedef unspecified_type Has_on_boundary_3;

  /*!
    a model of `Kernel::HasOnPositiveSide_3`
  */
  typedef unspecified_type Has_on_positive_side_3;

  /*!
    a model of `Kernel::HasOnNegativeSide_3`
  */
  typedef unspecified_type Has_on_negative_side_3;

  /*!
    a model of `Kernel::OrientedSide_3`
  */
  typedef unspecified_type Oriented_side_3;

  /*!
    a model of `Kernel::BoundedSide_3`
  */
  typedef unspecified_type Bounded_side_3;

  /*!
    a model of `Kernel::AreParallel_3`
  */
  typedef unspecified_type Are_parallel_3 ;

  /*!
    a model of `Kernel::AreOrderedAlongLine_3`
  */
  typedef unspecified_type Are_ordered_along_line_3 ;

  /*!
    a model of `Kernel::AreStrictlyOrderedAlongLine_3`
  */
  typedef unspecified_type Are_strictly_ordered_along_line_3;

  /*!
    a model of `Kernel::CollinearAreOrderedAlongLine_3`
  */
  typedef unspecified_type Collinear_are_ordered_along_line_3;

  /*!
    a model of `Kernel::CollinearAreStrictlyOrderedAlongLine_3`
  */
  typedef unspecified_type Collinear_are_strictly_ordered_along_line_3;

  /*!
    a model of `Kernel::DoIntersect_3`
  */
  typedef unspecified_type Do_intersect_3;

  /// @}

  /// \name Operations
  /// For each of the function objects above, there must exist a member
  /// function that requires no arguments and returns an instance of
  /// that function object. The name of the member function is the
  /// uncapitalized name of the type returned with the suffix `_object`
  /// appended. For example, for the function object
  /// `Kernel::Construct_vector_2` the following member function must
  /// exist:

  /// @{

  /*!

   */
  Kernel::Construct_vector_2 construct_vector_2_object() const ;

  /// @}

}; /* end Kernel */
