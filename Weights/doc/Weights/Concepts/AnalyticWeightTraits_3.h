/*!
\ingroup PkgWeightsRefConcepts
\cgalConcept

A concept that describes the set of requirements of classes used in the computation
of analytic weights in 3D.

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModelsEnd
*/
class AnalyticWeightTraits_3 {

public:

/// \name Types
/// @{

/*!
  A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/// @}

/// \name Geometric Objects
/// @{

/*!
  3D point type.
*/
typedef unspecified_type Point_3;

/*!
  3D vector type.
*/
typedef unspecified_type Vector_3;

/// @}

/// \name Constructions
/// @{

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Point_3& p, const Point_3& q)`

  that returns the squared Euclidean distance between the points `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_3;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Vector_3& v)`

  that returns the squared length of the vector `v`.
*/
typedef unspecified_type Compute_squared_length_3;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Vector_3& v, const Vector_3& w)`

  that returns the scalar product of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
  A construction object that must provide the function operator:

  `%Vector_3 operator()(const Vector_3& v, const Vector_3& w)`

  that returns the cross product of the vectors `v` and `w`.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/*!
  A construction object that must provide the function operator:

  `%Point_3 operator()(const Point_3& p, const Point_3& q, const Point_3& r)`

  that returns the center of the circle passing through the points `p`, `q`, and `r`.

  \pre The points `p`, `q`, and `r` are not collinear.
*/
typedef unspecified_type Construct_circumcenter_3;

/*!
  A construction object that must provide the function operator:

  `%Vector_3 operator()(const Point_3& p, const Point_3& q)`

  that returns the vector through the points `p` and `q`.
*/
typedef unspecified_type Construct_vector_3;

/*!
  A construction object that must provide the function operator:

  `%Point_3 operator()(const Point_3& p, const Point_3& q)`

  that returns the midpoint between the points `p` and `q`.
*/
typedef unspecified_type Construct_midpoint_3;

/*!
  A construction object that must provide two function operators:

  `%Point_3 operator()(const Point_3& p, const Point_3& q, const Point_3& r)`

  that returns the centroid of the points `p`, `q`, and `r` and

  `%Point_3 operator()(const Point_3& p, const Point_3& q, const Point_3& r, const Point_3& s)`

  that returns the centroid of the points `p`, `q`, `r` and `s`.
*/
typedef unspecified_type Construct_centroid_3;

/// @}

};
