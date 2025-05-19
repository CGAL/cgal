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

  /*! 2D point type */
  typedef unspecified_type Point_2;
  /*! 2D segment type */
  typedef unspecified_type Segment_2;
  /*! 2D triangle type */
  typedef unspecified_type Triangle_2;
  /*! 2D line type */
  typedef unspecified_type Line_2;
  /*! 2D vector type */
  typedef unspecified_type Vector_2;
  /*! 2D weighted point type */
  typedef unspecified_type Weighted_point_2;

  /*! 3D point type */
  typedef unspecified_type Point_3;
  /*! 3D line type */
  typedef unspecified_type Line_3;
  /*! 3D segment type */
  typedef unspecified_type Segment_3;
  /*! 3D vector type */
  typedef unspecified_type Vector_3;
  /*! 3D weighted point type */
  typedef unspecified_type Weighted_point_3;
  /*! 3D triangle type */
  typedef unspecified_type Triangle_3;
  /*! 3D ray type */
  typedef unspecified_type Ray_3;

  /*! `CGAL::Bounded_side` or `Uncertain<CGAL::Bounded_side>` */
  typedef unspecified_type Bounded_side;
  /*! `CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>` */
  typedef unspecified_type Comparison_result;
  /*! `bool` or `Uncertain<bool>` */
  typedef unspecified_type Boolean;

  /*! A function object to construct a `Point_2` */
  typedef unspecified_type Construct_point_2;
  /*! A function object to construct a `Weighted_point_2` */
  typedef unspecified_type Construct_weighted_point_2;
  /*! A function object to construct a `Segment_2` */
  typedef unspecified_type Construct_segment_2;
  /*! A function object to construct a `Vector_2` */
  typedef unspecified_type Construct_vector_2;
  /*! A function object to construct a `Triangle_2` */
  typedef unspecified_type Construct_triangle_2;
  /*! A function object to construct a `Line_2` */
  typedef unspecified_type Construct_line_2;

  /*! A function object to construct a translated `Point_2` */
  typedef unspecified_type Construct_translated_point_2;
  /*! A function object to construct the midpoint of two `Point_2` */
  typedef unspecified_type Construct_midpoint_2;
  /*! A function object to construct the barycenter of three `Point_2` */
  typedef unspecified_type Construct_barycenter_2;
  /*! A function object to construct the scaled version of a `Vector_2` */
  typedef unspecified_type Construct_scaled_vector_2;

  /*! A function object to construct the `i`-th vertex of a `Triangle_3` */
  typedef unspecified_type Construct_vertex_3;
  /*! A function object to construct the source vertex of a `Segment_3` */
  typedef unspecified_type Construct_source_3;
  /*! A function object to construct the target vertex of a `Segment_3` */
  typedef unspecified_type Construct_target_3;

  /*! A function object to compare signed distances to a line */
  typedef unspecified_type Compare_signed_distance_to_line_2;
  /*! A function object to compare signed distances to a line */
  typedef unspecified_type Less_signed_distance_to_line_2;

  /*! A function object to compute orientated side of a segment with respect to a triangle */
  typedef unspecified_type Oriented_side_2;
  /*! A function object to compare coordinates of two points in the `xy` order */
  typedef unspecified_type Compare_xy_2;
  /*! A function object to compare coordinates of two points from their `x` coordinates */
  typedef unspecified_type Less_x_3;
  /*! A function object to compare coordinates of two points from their `y` coordinates */
  typedef unspecified_type Less_y_3;
  /*! A function object to compare coordinates of two points from their `z` coordinates */
  typedef unspecified_type Less_z_3;
  /*! A function object to compare coordinates of two points from their `x` coordinates */
  typedef unspecified_type Compare_x_3;
  /*! A function object to compare coordinates of two points from their `y` coordinates */
  typedef unspecified_type Compare_y_3;
  /*! A function object to compare coordinates of two points from their `z` coordinates */
  typedef unspecified_type Compare_z_3;
  /*! A function object to compare coordinates of two points from their `x` coordinates */
  typedef unspecified_type Equal_x_3;
  /*! A function object to compare coordinates of two points from their `y` coordinates */
  typedef unspecified_type Equal_y_3;
  /*! A function object to compare coordinates of two points from their `y` coordinates */
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

/// @}

}; /* end IntersectionGeometricTraits_3 */
