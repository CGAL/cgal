/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

The concept `IsosurfacingTraits_3` describes the set of requirements to be
fulfilled by the traits class of a model of `IsosurfacingDomain_3`.

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the concept `Kernel`}
\cgalHasModelsEnd
*/
class IsosurfacingTraits_3
{
public:
  /// \name Types
  /// @{

  /*!
  The scalar type.
  Must be a model of `FieldNumberType`.
  */
  typedef unspecified_type FT;

  /*!
  The 3D point type.
  Must be a model of `Kernel::Point_3`.
  */
  typedef unspecified_type Point_3;

  /*!
  The 3D vector type.
  Must be a model of `Kernel::Vector_3`.
  */
  typedef unspecified_type Vector_3;

  /*!
  The 3D cuboid type.
  Must be a model of `Kernel::IsoCuboid_3`
  */
  typedef unspecified_type IsoCuboid_3;

  /*!
  A construction object that must provide the function operators:

  `FT operator()(Point_3 p)`

  and

  `FT operator()(Vector_3 p)`

  which return the \f$ x\f$-coordinate of the point and the vector, respectively.
  */
  typedef unspecified_type Compute_x_3;

  /*!
  A construction object that must provide the function operators:

  `FT operator()(Point_3 p)`

  and

  `FT operator()(Vector_3 p)`

  which return the \f$ y\f$-coordinate of the point and the vector, respectively.
  */
  typedef unspecified_type Compute_y_3;

  /*!
  A construction object that must provide the function operators:

  `FT operator()(Point_3 p)`

  and

  `FT operator()(Vector_3 p)`

  which return the \f$ z\f$-coordinate of the point and the vector, respectively.
  */
  typedef unspecified_type Compute_z_3;

  /*!
  A construction object that must provide the function operator:

  `Point_3 operator()(FT x, FT y, FT z)`

  which constructs a 3D point from its three coordinates.
  */
  typedef unspecified_type Construct_point_3;

  /*!
  A construction object that must provide the function operator:

  `Vector_3 operator()(FT x, FT y, FT z)`

  which constructs a 3D vector from its three coordinates.
  */
  typedef unspecified_type Construct_vector_3;

  /*!
  A construction object that must provide the function operator:

  `IsoCuboid_3 operator()(Point_3 p, Point_3 q)`

  which constructs an iso-oriented cuboid with diagonal opposite vertices `p` and `q`
  such that `p` is the lexicographically smallest point in the cuboid.
  */
  typedef unspecified_type Construct_iso_cuboid_3;

  /*!
  A construction object that must provide the function operator:

  `Point_3 operator()(IsoCuboid_3 c, int i)`

  which returns the i-th vertex of an iso-cuboid `c`. See `Kernel::ConstructVertex_3`
  for the order of vertices.
  */
  typedef unspecified_type Construct_vertex_3;

  /// @}

  /// \name Operations
  /// The following functions give access to the predicate and construction objects:
  /// @{

  /*!

  */
  Compute_x_3 compute_x_3_object();

  /*!

  */
  Compute_y_3 compute_y_3_object();

  /*!

  */
  Compute_z_3 compute_z_3_object();

  /*!

  */
  Construct_point_3 construct_point_3_object();

  /*!

  */
  Construct_vector_3 construct_vector_3_object();

  /*!

  */
  Construct_iso_cuboid_3 construct_iso_cuboid_3_object();

  /*!

  */
  Construct_vertex_3 construct_vertex_3_object();

  /// @}

};
