
/*!
\ingroup PkgInterpolation2Concepts
\cgalConcept

Most interpolation functions are parameterized by a traits class that
defines the primitives used in the interpolation algorithms. The concept
`InterpolationTraits` defines this common set of requirements.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Interpolation_traits_2<K>}
\cgalHasModels{CGAL::Interpolation_gradient_fitting_traits_2<K>}
\cgalHasModelsEnd

\sa `GradientFittingTraits`
\sa `CGAL::sibson_c1_interpolation()`
\sa \ref PkgInterpolationSibsonGradientFitting
\sa `CGAL::farin_c1_interpolation()`
\sa `CGAL::quadratic_interpolation()`

*/
class InterpolationTraits
{
public:

/// \name Types
/// @{

/*!
The number type must follow the model `FieldNumberType`.
*/
typedef unspecified_type FT;

/*!
The (weightless) point type.
*/
typedef unspecified_type Point_d;

/*!
The weighted point type.
*/
typedef unspecified_type Weighted_point_d;

/*!
The corresponding vector type.
*/
typedef unspecified_type Vector_d;

/*!
A constructor object for `Point_d`. Provides :

`Point_d operator() (Point_d p)` which returns the point itself.

`Point_d operator() (Weighted_point_d wp)` which extracts the bare point from the weighted point.
*/
typedef unspecified_type Construct_point_d;

/*!
A constructor object for `Vector_d`. Provides :

`Vector_d operator() (Point_d p, Point_d q)` which produces the
vector `q - p` and

`Vector_d operator() (Null_vector NULL_VECTOR)` which introduces
the null vector.
*/
typedef unspecified_type Construct_vector_d;

/*!
Constructor object for `Vector_d`. Provides :

`Vector_d operator() (Vector_d v, FT scale)` which produces the
vector `v` scaled by a factor `scale`.
*/
typedef unspecified_type Construct_scaled_vector_d;

/*!
Constructor object for `FT`. Provides the operator:

`FT operator() (Point_d p, Point_d q)` returning the squared
distance between `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_d;

/*!
default constructor.
*/
InterpolationTraits();

/// @}

/// \name Construction objects
/// The following functions that create instances of the above
/// constructor object types must exist.
/// @{

/*!

*/
Construct_vector_d construct_point_d_object();

/*!

*/
Construct_vector_d construct_vector_d_object();

/*!

*/
Construct_scaled_vector_d construct_scaled_vector_d_object();

/*!

*/
Compute_squared_distance_d compute_squared_distance_d_object();

/// @}

}; /* end InterpolationTraits */

