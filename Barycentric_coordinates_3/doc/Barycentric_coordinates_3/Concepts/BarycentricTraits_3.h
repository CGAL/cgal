namespace CGAL {
namespace Barycentric_coordinates {

/*!
\ingroup PkgBarycentricCoordinates3RefConcepts
\cgalConcept

A concept that describes the set of requirements of the template parameter
`GeomTraits` used to parameterize all classes and functions with 3D barycentric
coordinates from the namespace `CGAL::Barycentric_coordinates`. The concept `BarycentricTraits_3`
extends the concept `BarycentricTraits_2` and adds the requirements for 3D objects.
\cgalGeneralizes `BarycentricTraits_2`

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModelsEnd

*/
class BarycentricTraits_3 {

public:

/// \name Types
/// @{

/*!
  A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/// @}

/// \name 3D Geometric Objects
/// @{

/*!
  A model of `Kernel::Point_3`.
*/
typedef unspecified_type Point_3;

/*!
  A model of `Kernel::Vector_3`.
*/
typedef unspecified_type Vector_3;

/*!
  A model of `Kernel::Plane_3`.
*/
typedef unspecified_type Plane_3;

/// @}

/// \name 2D Geometric Objects
/// @{

/*!
  A model of `Kernel::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
  A model of `Kernel::Vector_2`.
*/
typedef unspecified_type Vector_2;

/*!
  A model of `Kernel::Triangle_2`.
*/
typedef unspecified_type Triangle_2;

/*!
  A model of `Kernel::Iso_rectangle_2`.
*/
typedef unspecified_type Iso_rectangle_2;

/// @}

/// \name 3D Generalized Constructions
/// @{

/*!
  A construction object that must provide the function operator:

  `FT operator(Vector_3 u, Vector_3 v)`

  that returns an approximation of the angle between `u` and `v`.
  The angle is given in degrees.
*/
typedef unspecified_type Compute_approximate_angle_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(Vector_3 u, Vector_3 v, Vector_3 w)`

  that returns the determinant of the three vectors `u`, `v` and `w`.
*/
typedef unspecified_type Compute_determinant_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(Vector_3 v, Vector_3 w)`

  that returns the scalar (inner) product of the two vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
  A construction object that must provide the function operator:

  `FT operator()(const Vector_3& v)`

  that returns the squared length of the vector `v`.
*/
typedef unspecified_type Compute_squared_length_3;

/*!
  A construction object that must provide the function operator:

  `FT operator(Point_3 p0, Point_3 p1, Point_3 p2, Point_3 p3)`

  that returns the signed volume of the tetrahedron defined by the four points `p0`, `p1`, `p2`, and `p3`.
*/
typedef unspecified_type Compute_volume_3;

/*!
  A construction object that must provide the function operator:

  `Vector_3 operator(Vector_3 u, Vector_3 v)`

  that returns the cross product between `u` and `v`.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/*!
  A construction object that must provide the function operator:

  `Vector_3 operator()(Vector_3 v, FT t) const`

  returning the vector `v / t`.
*/
typedef unspecified_type Construct_divided_vector_3;

/*!
  A construction object that must provide the function operator:

  `Vector_3 operator(Point_3 p, Point_3 q)`

  that returns the vector `q` - `p`.
*/
typedef unspecified_type Construct_vector_3;
/// @}

/// \name 2D Generalized Constructions
/// @{

/*!
  A construction object that must provide the function operator:

  `FT operator(Point_2 p0, Point_2 p1, Point_3 p2)`

  that returns the signed area of the triangle defined by the three points `p0`, `p1`, and `p2`.
*/
typedef unspecified_type Compute_area_2;

/*!
  A construction object that must provide the function operator:

 `FT operator()(Point_2 p, Point_2 q),`

  which returns the squared distance between `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_2;

/*!
  A construction object that must provide the function operator:

  `FT operator()(Vector_2 v, Vector_2 w)`

  that returns the scalar product of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_2;

/*!
  A construction object that must provide the function operator:

  `Vector_2 operator()(Point_2 a, Point_2 b)`

  that computes the vector `b-a`.
*/
typedef unspecified_type Construct_vector_2;
/// @}

/// \name 2D Generalized Predicates
/// @{

/*!
  A predicate object that must provide the function operator:

  `bool operator()(Point_2 p, Point_2 q, Point_2 r)`

  that returns `true` if the points `p`, `q`, and `r` are collinear and `false` otherwise.
*/
typedef unspecified_type Collinear_2;

/*!
  A predicate object that must provide the function operator:

  `bool operator()(Point_2 p, Point_2 q, Point_2 r)`

  that returns `true`, iff `q` lies between `p` and `r` and `p`, `q`,
  and `r` satisfy the precondition that they are collinear.
*/
typedef unspecified_type Collinear_are_ordered_along_line_2;

/*!
  A predicate object that must provide the function operator:

  `Comparison_result operator()(Point_2 p, Point_2 q)`

  that returns `SMALLER, EQUAL` or `LARGER`
  according to the \f$ x\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
  A predicate object that must provide the function operator:

  `Comparison_result operator()(Point_2 p, Point_2 q)`

  that returns `SMALLER, EQUAL` or `LARGER`
  according to the \f$ y\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;

/*!
  A predicate object that must provide the function operator:

  `bool operator()(Point_2 p, Point_2 q)`,

  which returns `true` iff `p` is equal to `q`.
*/
typedef unspecified_type Equal_2;

/*!
  A predicate object that must provide the function operator:

  `bool operator()(Point_2 p, Point_2 q)`

  where `true` is returned iff \f$ p <_{xy} q\f$.
  We have \f$ p<_{xy}q\f$, iff \f$ p_x < q_x\f$ or \f$ p_x = q_x\f$ and \f$ p_y < q_y\f$,
  where \f$ p_x\f$ and \f$ p_y\f$ denote \f$ x\f$ and \f$ y\f$ coordinate of point \f$ p\f$,
  respectively.
*/
typedef unspecified_type Less_xy_2;

/*!
  A predicate object that must provide the function operator:

  `Orientation operator()(Point_2 p, Point_2 q, Point_2 r)`

  that returns `CGAL::LEFT_TURN` if `r` lies to the left of the oriented line `l` defined by `p` and `q`,
  returns `CGAL::RIGHT_TURN` if `r` lies to the right of `l`, and returns `CGAL::COLLINEAR` if `r` lies on `l`.
*/
typedef unspecified_type Orientation_2;
/// @}

/// The following functions give access to the predicate and constructor objects.
/// @{

/*!

*/
Collinear_2 collinear_2_object();

/*!

*/
Collinear_are_ordered_along_line_2 collinear_are_ordered_along_line_2_object();

/*!

*/
Compute_approximate_angle_3 compute_approximate_angle_3_object();

/*!

*/
Compute_area_2 compute_area_2_object();

/*!

*/
Compute_determinant_3 compute_determinant_3_object();

/*!

*/
Compute_scalar_product_2 compute_scalar_product_2_object();

/*!

*/
Compute_scalar_product_3 compute_scalar_product_3_object();

/*!

*/
Compute_squared_distance_2 compute_squared_distance_2_object();

/*!

*/
Compute_volume_3 compute_volume_3_object();

/*!

*/
Compare_x_2 compare_x_2_object();

/*!

*/
Compare_y_2 compare_y_2_object();

/*!

*/
Construct_cross_product_vector_3 construct_cross_product_vector_3_object();

/*!

*/
Construct_divided_vector_3 construct_divided_vector_3_object();

/*!

*/
Construct_vector_2 construct_vector_2_object();

/*!

*/
Construct_vector_3 construct_vector_3_object();

/*!

*/
Equal_2 equal_2_object();

/*!

*/
Less_xy_2 less_xy_2_object();

/*!

*/
Orientation_2 orientation_2_object();
/// @}

};

} // namespace Barycentric_coordinates
} // namespace CGAL
