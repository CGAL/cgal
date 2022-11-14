/*!
\ingroup PkgWeightsRefConcepts
\cgalConcept

A concept that describes the set of requirements of classes used in the computation
of analytic weights in 2D.

\cgalHasModel
- All models of `Kernel`
- `CGAL::Projection_traits_xy_3<K>`
- `CGAL::Projection_traits_yz_3<K>`
- `CGAL::Projection_traits_xz_3<K>`
*/
class AnalyticWeightTraits_2 {

public:

/// \name Types
/// @{

/*!
  A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/*!
  `CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>`.
*/
typedef unspecified_type Comparison_result;

/*!
  `CGAL::Orientation` or `Uncertain<CGAL::Orientation>`.
*/
typedef unspecified_type Orientation;

/// @}

/// \name Geometric Objects
/// @{

/*!
  2D point type.
*/
typedef unspecified_type Point_2;

/*!
  2D vector type.
*/
typedef unspecified_type Vector_2;

/// @}

/// \name Constructions
/// @{

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Point_2& p, const Point_2& q, const Point_2& r)`

  that returns the signed area of the triangle defined by the points `p`, `q`, and `r`.
*/
typedef unspecified_type Compute_area_2;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Point_2& p, const Point_2& q)`

  that returns the squared Euclidean distance between the points `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_2;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Vector_2& v)`

  that returns the squared length of the vector `v`.
*/
typedef unspecified_type Compute_squared_length_2;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Vector_2& v, const Vector_2& w)`

  that returns the scalar product of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_2;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Vector_2& v, const Vector_2& w)`

  that returns the determinant of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_determinant_2;

/*!
  A construction object that must provide the function operator:

  `%Point_2 operator()(const Point_2& p, const Point_2& q, const Point_2& r)`

  that returns the center of the circle passing through the points `p`, `q`, and `r`.

  \pre The points `p`, `q`, and `r` are not collinear.
*/
typedef unspecified_type Construct_circumcenter_2;

/*!
  A construction object that must provide the function operator:

  `%Vector_2 operator()(const Point_2& p, const Point_2& q)`

  that returns the vector through the points `p` and `q`.
*/
typedef unspecified_type Construct_vector_2;

/*!
  A construction object that must provide the function operator:

  `%Point_2 operator()(const Point_2& p, const Point_2& q)`

  that returns the midpoint between the points `p` and `q`.
*/
typedef unspecified_type Construct_midpoint_2;

/*!
  A construction object that must provide the function operator:

  `%Point_2 operator()(const Point_2& p, const Point_2& q, const Point_2& r)`

  that returns the centroid of the points `p`, `q`, and `r`.
*/
typedef unspecified_type Construct_centroid_2;

/// @}

/// \name Predicates
/// @{

/*!
  A predicate object that must provide the function operator:

  `bool operator()(const Point_2& p, const Point_2& q)`

  that returns `true` if `p = q` and `false` otherwise.
*/
typedef unspecified_type Equal_2;

/*!
  A predicate object that must provide the function operator:

  `bool operator()(const Point_2& p, const Point_2& q, const Point_2& r)`

  that returns `true` if the points `p`, `q`, and `r` are collinear and `false` otherwise.
*/
typedef unspecified_type Collinear_2;

/*!
  A predicate object that must provide the function operator:

  `bool operator()(const Point_2& p, const Point_2& q)`

  that returns `true` iff the x-coordinate of `p` is smaller than the x-coordinate of `q` or
  if they are the same and the y-coordinate of `p` is smaller than the y-coordinate of `q`.
*/
typedef unspecified_type Less_xy_2;

/*!
  A predicate object that must provide the function operator:

  `Comparison_result operator()(const Point_2& p, const Point_2& q)`

  that compares the %Cartesian x-coordinates of the points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
  A predicate object that must provide the function operator:

  `Comparison_result operator()(const Point_2& p, const Point_2& q)`

  that compares the %Cartesian y-coordinates of the points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;

/*!
  A predicate object that must provide the function operator:

  `Orientation operator()(const Point_2& p, const Point_2& q, const Point_2& r)`

  that returns `CGAL::LEFT_TURN` if `r` lies to the left of the oriented line `l` defined by `p` and `q`,
  returns `CGAL::RIGHT_TURN` if `r` lies to the right of `l`, and returns `CGAL::COLLINEAR` if `r` lies on `l`.
*/
typedef unspecified_type Orientation_2;

/// @}

};
