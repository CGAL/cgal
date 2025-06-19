/*!
\cgalConcept
\ingroup PkgKernel23Concepts

The concept `ProjectionTraitsGeometricTraits_3` provides types and functors
required to be used with the class `CGAL::Projection_traits_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::Projection_traits_3`

*/
class ProjectionTraitsGeometricTraits_3
{
public:

/// \name Types
/// @{

  /*! Number type model of `FieldNumberType`*/
  typedef unspecified_type FT;

  /*! 2D line type */
  typedef unspecified_type Line_2;
  /*! 2D point type */
  typedef unspecified_type Point_2;
  /*! 2D segment type */
  typedef unspecified_type Segment_2;
  /*! 2D triangle type */
  typedef unspecified_type Triangle_2;
  /*! 2D vector type */
  typedef unspecified_type Vector_2;
  /*! 2D weighted point type */
  typedef unspecified_type Weighted_point_2;

  /*! 3D point type */
  typedef unspecified_type Point_3;
  /*! 3D line type */
  typedef unspecified_type Line_3;
  /*! 3D ray type */
  typedef unspecified_type Ray_3;
  /*! 3D segment type */
  typedef unspecified_type Segment_3;
  /*! 3D triangle type */
  typedef unspecified_type Triangle_3;
  /*! 3D vector type */
  typedef unspecified_type Vector_3;
  /*! 3D weighted point type */
  typedef unspecified_type Weighted_point_3;

  /*! `CGAL::Bounded_side` or `Uncertain<CGAL::Bounded_side>` */
  typedef unspecified_type Bounded_side;
  /*! `CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>` */
  typedef unspecified_type Comparison_result;
  /*! `bool` or `Uncertain<bool>` */
  typedef unspecified_type Boolean;

  /*! A functor model of `ConstructPoint_3` */
  typedef unspecified_type Construct_point_3;
  /*! A functor model of `ConstructWeightedPoint_3` */
  typedef unspecified_type Construct_weighted_point_3;
  /*! A functor model of `ConstructSegment_3` */
  typedef unspecified_type Construct_segment_3;
  /*! A functor model of `ConstructVector_3` */
  typedef unspecified_type Construct_vector_3;
  /*! A functor model of `ConstructTriangle_3` */
  typedef unspecified_type Construct_triangle_3;
  /*! A functor model of `ConstructLine_3` */
  typedef unspecified_type Construct_line_3;
  /*! A functor model of `ConstructTriangle_3` */
  typedef unspecified_type Construct_triangle_3;
  /*! A functor model of `ConstructBbox_3` */
  typedef unspecified_type Construct_bbox_3;

  /*! A functor model of `ConstructTranslatedPoint_3` */
  typedef unspecified_type Construct_translated_point_3;
  /*! A functor model of `ConstructMidpoint_3` */
  typedef unspecified_type Construct_midpoint_3;
  /*! A functor model of `ConstructBarycenter_3` */
  typedef unspecified_type Construct_barycenter_3;
  /*! A functor model of `ConstructVector_3` */
  typedef unspecified_type Construct_vector_3;
  /*! A functor model of `ConstructScaledVector_3` */
  typedef unspecified_type Construct_scaled_vector_3;

  /*! A functor model of `ConstructVertex_3` */
  typedef unspecified_type Construct_vertex_3;
  /*! A functor model of `ConstructSource_3` */
  typedef unspecified_type Construct_source_3;
  /*! A functor model of `ConstructTarget_3` */
  typedef unspecified_type Construct_target_3;

  /*! A functor model of `CompareSignedDistanceToLine_2` */
  typedef unspecified_type Compare_signed_distance_to_line_2;
  /*! A functor model of `LessSignedDistanceToLine_2` */
  typedef unspecified_type Less_signed_distance_to_line_2;

  /*! A functor model of `OrientedSide_2` */
  typedef unspecified_type Oriented_side_2;
  /*! A functor model of `CompareXYZ_3` */
  typedef unspecified_type Compare_xyz_3;
  /*! A functor model of `LessX_3` */
  typedef unspecified_type Less_x_3;
  /*! A functor model of `LessY_3` */
  typedef unspecified_type Less_y_3;
  /*! A functor model of `LessZ_3` */
  typedef unspecified_type Less_z_3;
  /*! A functor model of `CompareX_3` */
  typedef unspecified_type Compare_x_3;
  /*! A functor model of `CompareY_3` */
  typedef unspecified_type Compare_y_3;
  /*! A functor model of `CompareZ_3` */
  typedef unspecified_type Compare_z_3;
  /*! A functor model of `EqualX_3` */
  typedef unspecified_type Equal_x_3;
  /*! A functor model of `EqualY_3` */
  typedef unspecified_type Equal_y_3;
  /*! A functor model of `EqualZ_3` */
  typedef unspecified_type Equal_z_3;

/// @}

/// \name Operations
/// @{

/*!
* returns a function object model of `ComputeArea_2`
*/
unspecified_type compute_area_2_object();

/*!
* returns a function object model of `ConstructPoint_3`
*/
unspecified_type construct_point_3_object();

/*!
 * returns a function object model of `Compute_SquaredRadius_2`
 */
unspecified_type compute_squared_radius_2_object();

/// @}

}; /* end IntersectionGeometricTraits_3 */
