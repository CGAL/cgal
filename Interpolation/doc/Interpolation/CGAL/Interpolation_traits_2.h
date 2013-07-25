
namespace CGAL {

/*!
\ingroup PkgInterpolation2Interpolation

`Interpolation_traits_2` is a model for the concept 
`InterpolationTraits` and can be used to instantiate the 
geometric traits class of interpolation methods applied on a 
bivariate function over a two-dimensional domain. The traits class 
is templated by a kernel class `K`. 

\cgalModels `InterpolationTraits`

\sa `InterpolationTraits` 
\sa `GradientFittingTraits` 
\sa `CGAL::Interpolation_gradient_fitting_traits_2<K>` 

*/
template< typename K >
class Interpolation_traits_2 {
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
typedef K::Construct_vector_2 Construct_vector_d; 

/*!

*/ 
typedef K::Construct_scaled_vector_2 
Construct_scaled_vector_d; 

/*!

*/ 
typedef K::Compute_squared_distance_2 
Compute_squared_distance_d; 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
Construct_scaled_vector_d 
construct_scaled_vector_d_object() const; 

/*!

*/ 
Construct_vector_d 
construct_vector_d_object()const; 

/*!

*/ 
Compute_squared_distance_d 
compute_squared_distance_d_object() const; 

/// @}

}; /* end Interpolation_traits_2 */
} /* end namespace CGAL */
