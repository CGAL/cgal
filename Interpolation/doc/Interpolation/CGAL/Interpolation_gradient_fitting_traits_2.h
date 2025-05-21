
namespace CGAL {

/*!
\ingroup PkgInterpolation2Interpolation

`Interpolation_gradient_fitting_traits_2` is a model for the concepts
`InterpolationTraits` and `GradientFittingTraits`. It can be
used to instantiate the geometric traits class of interpolation
functions and of Sibson's gradient fitting function when applied on
a function defined over a two-dimensional domain. The traits class
is templated by a kernel class `K`.

\cgalModels{GradientFittingTraits,InterpolationTraits}

\sa `InterpolationTraits`
\sa `GradientFittingTraits`
\sa `CGAL::Interpolation_traits_2<K>`

*/
template< typename K >
class Interpolation_gradient_fitting_traits_2
{
public:

/// \name Types
/// @{

/*!

*/
typedef K::FT FT;

/*!

*/
typedef K::Point_2 Point_d;

/*!

*/
typedef K::Vector_2 Vector_d;

/*!

*/
typedef K::Aff_transformation_2 Aff_transformation_d;

/*!

*/
typedef K::Construct_vector_2 Construct_vector_d;

/*!

*/
typedef K::Construct_scaled_vector_2 Construct_scaled_vector_d;

/*!

*/
typedef K::Compute_squared_distance_2 Compute_squared_distance_d;

/*!

*/
typedef Construct_null_matrix_2<Aff_transformation_d> Construct_null_matrix_d;

/*!

*/
typedef Construct_scaling_matrix_2<Aff_transformation_d> Construct_scaling_matrix_d;

/*!

*/
typedef Construct_sum_matrix_2<Aff_transformation_d> Construct_sum_matrix_d;

/*!

*/
typedef Construct_outer_product_2<K> Construct_outer_product_d;

/// @}

/// \name Operations
/// @{

/*!

*/
Construct_scaled_vector_d construct_scaled_vector_d_object() const;

/*!

*/
Construct_vector_d construct_vector_d_object()const;

/*!

*/
Compute_squared_distance_d compute_squared_distance_d_object() const;

/*!

*/
Construct_null_matrix_d construct_null_matrix_d_object() const;

/*!

*/
Construct_scaling_matrix_d construct_scaling_matrix_d_object() const;

/*!

*/
Construct_sum_matrix_d construct_sum_matrix_d_object() const;

/*!

*/
Construct_outer_product_d construct_outer_product_d_object() const;

/// @}

}; /* end Interpolation_gradient_fitting_traits_2 */
} /* end namespace CGAL */
