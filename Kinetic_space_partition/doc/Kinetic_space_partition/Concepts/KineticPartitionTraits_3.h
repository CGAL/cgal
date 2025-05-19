 /*!
\ingroup PkgKineticSpacePartitionConcepts
\cgalConcept

A concept that describes the set of types required by the `CGAL::Kinetic_space_partition_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd

\sa `CGAL::Kinetic_space_partition_3`
*/
class KineticSpacePartitionTraits_3 {

public:

/// \name Types
/// @{

  /// The 3D point type
  typedef unspecified_type Point_3;
  /// The 2D point type
  typedef unspecified_type Point_2;
  /// The 3D vector type
  typedef unspecified_type Vector_3;
  /// The 2D line type
  typedef unspecified_type Line_2;
  /// The 2D direction type
  typedef unspecified_type Direction_2;
  /// The plane type
  typedef unspecified_type Plane_3;
  /// The 2D vector type
  typedef unspecified_type Vector_2;
  /// The 2D segment type
  typedef unspecified_type Segment_2;
  /// The 3d tetrahedron type
  typedef unspecified_type Tetrahedron_3;
  /// The 3d transformation type
  typedef unspecified_type Transform_3;

  /// The number type of the Cartesian coordinates
  typedef unspecified_type FT;

  /*!
   * Function object type that provides
   * `Point_3 operator()(Origin p)`
   * returning the point with 0, 0, 0 as Cartesian coordinates
   * and `Point_3 operator()(FT x, FT y, FT z)`
   * returning the point with `x`, `y` and `z` as Cartesian coordinates.
   */
  typedef unspecified_type Construct_point_3;

  /*!
   * Function object type that provides
   * `Vector_3 operator()(Point_3 p1, Point_3 p2)` and
   * `Vector_3 operator()(Origin p1, Point_3 p2)`
   * returning the vector `p1p2`, `Vector_3 operator()(NULL_VECTOR)` returning
   * the null vector, and `Vector_3 operator()(Line_3 l)` returning
   * a vector having the same direction as `l`
   */
  typedef unspecified_type Construct_vector_3;

  /*!
   * Function object type that provides `Line_2 operator()(Point_2 p, Vector_2 d)`
   * returning the line going through  `p` in the direction of `d`.
   */
  typedef unspecified_type Construct_line_2;

  /*!
   * Function object type that provides
   * `Point_2 operator()(Line_2 l, int i)`
   * returning an arbitrary point on `l`. `i` is not used and can be of
   * any value.
   */
  typedef unspecified_type Construct_point_on_2;

  /*!
   * Function object type that provides
   * `Point_2 operator()(FT x, FT y)`
   * returning the 2D point with `x` and `y` as Cartesian coordinates.
   */
  typedef unspecified_type Construct_point_2;

  /*!
   * Function object type that provides
   * `Vector_2 operator()(Point_2 p1, Point_2 p2)`
   * returning the vector `p1p2`, `Vector_2 operator()(NULL_VECTOR)` returning
   * the null vector.
   */
  typedef unspecified_type Construct_vector_2;

  /*!
   * Function object type that provides
   * `Tetrahedron_3 operator(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`
   * returning the tetrahedron with the points `p`, `q`, `r` and `s`.
   */
  typedef unspecified_type ConstructTetrahedron_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v)`
   * returning the `x` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_x_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v)`
   * returning the `y` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_y_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v)`
   * returning the `z` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_z_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_2 p)` and `FT operator()(Vector_2 v)`
   * returning the `x` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_x_2;

  /*!
   * Function object type that provides
   * `FT operator()(Point_2 p)` and `FT operator()(Vector_2 v)`
   * returning the `y` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_y_2;

  /*!
   * Function object type that provides
   * `FT operator()(Vector_2 v)`
   * returning the squared length of `v`.
   */
  typedef unspecified_type Compute_squared_length_2;

  /*!
   * Function object type that provides
   * `FT operator()(Vector_3 v1, Vector_3 v2)`
   * returning the scalar product of `v1` and `v2`.
   */
  typedef unspecified_type Compute_scalar_product_3;

  /*!
   * Function object type that provides
   * `Vector_3 operator() (Vector_3 v1, Vector_3 v2)`
   * returning the `v1+v2`.
   */
  typedef unspecified_type Construct_sum_of_vectors_3;

  /*!
   * Function object type that provides
   * `Vector_3 operator()(Plane_3 p)`
   * returns a vector that is orthogonal to the plane `p` and directed to the positive side of `p`.
   */
  typedef unspecified_type Construct_orthogonal_vector_3;

  /*!
   * Function object type that provides
   * `Plane_3 operator()(Point_3 p, Point_3 q, Point_3 r)`
   *  creates a plane passing through the points `p`, `q` and `r` and
   * `Plane_3 operator()(Point_3 p, Vector_3 v)`
   *  introduces a plane that passes through point `p` and that is orthogonal to `V`.
   */
  typedef unspecified_type Construct_plane_3;

  /*!
   * Function object type that provides
   * `Vector_3 operator()(Plane_3 h, Point_3 p)
   * returns the orthogonal projection of `p` onto `h`.
   */
  typedef unspecified_type Construct_projected_point_3;

  /*!
   * Function object type that provides
   * `bool operator()(Point_3 p, Point_3 q, Point_3 r)`
   * returning true if the points `p`, `q`, and `r` are collinear and
   * false otherwise.
   */
  typedef unspecified_type Collinear_3;

  /*!
   * Function object type that provides
   * `Oriented_size operator()(Plane_3 h, Point_3 p)`
   * returns `CGAL::ON_ORIENTED_BOUNDARY`, `CGAL::ON_NEGATIVE_SIDE`, or `CGAL::ON_POSITIVE_SIDE`, depending on the position of `p` relative to the oriented plane `h`.
   */
  typedef unspecified_type Oriented_side_3;

/// @}

/// \name Access to Function Objects
/// @{

  Construct_point_3
  construct_point_3_object();

  Construct_vector_3
  construct_vector_3_object();

  Construct_line_2
  construct_line_2_object();

  Construct_point_on_3
  construct_point_on_3_object();

  Construct_point_2
  construct_point_2_object();

  Construct_vector_2
  construct_vector_2_object();

  Construct_tetrahedron_3
  construct_tetrahedron_3_object();

  Compute_x_3
  compute_x_3_object();

  Compute_y_3
  compute_y_3_object();

  Compute_z_3
  compute_z_3_object();

  Compute_x_2
  compute_x_2_object();

  Compute_y_2
  compute_y_2_object();

  Compute_squared_length_2
  compute_squared_length_2_object();

  Construct_sum_of_vectors_3
  construct_sum_of_vectors_3_object();

  Construct_projected_point_3
  construct_projected_point_3_object();

  Compute_scalar_product_3
  compute_scalar_product_3_object();

  Collinear_3
  collinear_3_object();

  Oriented_side_3
  oriented_side_3_object();

/// @}
};
