
namespace CGAL {

/*!
\ingroup PkgInterpolation2

The struct `Data_access` implements a functor that allows to retrieve 
data from an associative container. The functor keeps a reference to 
the container. Given an instance of the container's key type, it 
returns a pair of the container's value type and a Boolean indicating 
whether the retrieval was successful. 

This class can be used to provide the values and gradients of the 
interpolation functions. 

Parameters 
-------------- 

The class 
`Data_access` has the container type `Map` as template parameter. 

*/
template< typename Map >
class Data_access {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Map::mapped_type Data_type; 

/*! 

*/ 
typedef Map::key_type Key_type; 

/// @} 

/// \name Creation 
/// @{

/*! 
Introduces a `Data_access` to the container `map`. 
*/ 
Data_access(const Map& map); 

/*! 
If 
there is an entry for `p` in the container `map`, then the 
pair of `map.find(p)` and `true` is returned. Otherwise, the 
Boolean value of the pair is `false`. 
*/ 
std::pair< Data_type, bool> operator()(const Key_type& p); 

/// @}

}; /* end Data_access */

/*!
\ingroup PkgInterpolation2Interpolation

generates the interpolated function value
computed by Farin's interpolant \cite f-sodt-90. See also
`sibson_c1_interpolation`. \pre `norm` \f$ \neq0\f$. `function_value(p).second == true` for all points `p` of the point/coordinate pairs in the range \f$ \left[\right.\f$`first`, `beyond`\f$ \left.\right)\f$.
\pre The range \f$ \left[\right.\f$ `first`, `beyond`\f$ \left.\right)\f$ contains either one or more than three element
The function `farin_c1_interpolation` interpolates the function values and the 
gradients that are provided by functors using the method described in \cite f-sodt-90. 

Parameters 
-------------- 

`RandomAccessIterator::value_type` is a pair 
associating a point to a (non-normalized) barycentric coordinate. See 
`sibson_c1_interpolation` for the other parameters. 

Requirements 
-------------- 

Same requirements as for `sibson_c1_interpolation` only the 
iterator must provide random access and `Traits::FT` does not need 
to provide the square root operation. 

\sa CGAL::Data_access<Map> 
\sa CGAL::linear_interpolation 
\sa CGAL::sibson_c1_interpolation 
\sa CGAL::sibson_gradient_fitting 
\sa CGAL::Interpolation_traits_2<K> 
\sa CGAL::natural_neighbor_coordinates_2 
\sa CGAL::regular_neighbor_coordinates_2 
\sa CGAL::surface_neighbor_coordinates_3 

s. 
*/
template < class RandomAccessIterator, class Functor,
class GradFunctor, class Traits> typename Functor::result_type
farin_c1_interpolation(RandomAccessIterator first,
RandomAccessIterator beyond, const typename
std::iterator_traits<RandomAccessIterator>::value_type::second_type&
norm, const typename
std::iterator_traits<ForwardIterator>::value_type::first_type& p,
Functor function_value, GradFunctor function_gradient, const
Traits& traits);


/*!
\ingroup PkgInterpolation2Interpolation

The function `linear_interpolation` computes the weighted sum of the function 
values which must be provided via a functor. 

`ForwardIterator::value_type` is a pair associating a point to a (non-normalized) barycentric
coordinate. `norm` is the normalization factor. Given a point,
the functor `function_values` allows to access a pair of a
function value and a Boolean. The Boolean indicates whether the
function value could be retrieved correctly. This function generates
the interpolated function value as the weighted sum of the values
corresponding to each point of the point/coordinate pairs in the
range \f$ \left[\right.\f$`first`, `beyond`\f$ \left.\right)\f$.
\pre `norm` \f$ \neq0\f$. `function_value(p).second == true` for all points `p` of the point/coordinate pairs in the range \f$ \left[\right.\f$`first`, `beyond`\f$ \left.\right)\f$.

Requirements 
-------------- 

<OL> 
<LI>`ForwardIterator::value_type` is a pair of 
point/coordinate value, thus 
`ForwardIterator::value_type::first_type` is equivalent to a 
point and `ForwardIterator::value_type::second_type` is a 
field number type. 
<LI>`Functor::argument_type` must be equivalent to 
`ForwardIterator::value_type::first_type` and 
`Functor::result_type` is a pair of the function value type 
and a Boolean value. The function value type must provide a 
multiplication and addition operation with the field number type 
`ForwardIterator::value_type::second_type` and a constructor 
with argument \f$ 0\f$. A model of the functor is provided by the 
struct `Data_access`. It must be instantiated accordingly with 
an associative container (e.g. \stl `std::map`) having the 
point type as `key_type` and the function value type as 
`mapped_type`. 
</OL> 

\sa CGAL::Data_access<Map> 
\sa CGAL::natural_neighbor_coordinates_2 
\sa CGAL::regular_neighbor_coordinates_2 
\sa CGAL::surface_neighbor_coordinates_3 

*/
template < class ForwardIterator, class Functor> typename
Functor::result_type::first_type linear_interpolation(ForwardIterator first,
ForwardIterator beyond, const typename
std::iterator_traits<ForwardIterator>::value_type::second_type&
norm, Functor function_values);


/*!
\ingroup PkgInterpolation2Interpolation

The function `quadratic_interpolation` interpolates the function values and first degree 
functions defined from the function gradients. Both, function values and 
gradients, must be provided by functors. 

This function generates the
interpolated function value as the weighted sum of the values plus a
linear term in the gradient for each point of the point/coordinate
pairs in the range \f$ \left[\right.\f$ `first`,
`beyond`\f$ \left.\right)\f$. See also
`sibson_c1_interpolation`. \pre `norm` \f$ \neq0\f$ `function_value(p).second == true` for all points `p` of the point/coordinate pairs in the range \f$ \left[\right.\f$`first`, `beyond`\f$ \left.\right)\f$.

Parameters 
-------------- 

See `sibson_c1_interpolation`. 

Requirements 
-------------- 

Same requirements as for 
`sibson_c1_interpolation` only that `Traits::FT` does not need 
to provide the square root operation. 

\sa `InterpolationTraits` 
\sa `GradientFittingTraits` 
\sa CGAL::Data_access<Map> 
\sa CGAL::sibson_gradient_fitting 
\sa CGAL::linear_interpolation 
\sa CGAL::Interpolation_traits_2<K> 
\sa CGAL::Interpolation_gradient_fitting_traits_2<K> 
\sa CGAL::natural_neighbor_coordinates_2 
\sa CGAL::regular_neighbor_coordinates_2 
\sa CGAL::surface_neighbor_coordinates_3 
*/
template < class ForwardIterator, class Functor, class
GradFunctor, class Traits> typename Functor::result_type
quadratic_interpolation(ForwardIterator first, ForwardIterator
beyond, const typename std::iterator_traits<ForwardIterator>::
value_type::second_type& norm, 
const typename std::iterator_traits<ForwardIterator>::value_type::
first_type& p, Functor function_value, GradFunctor
function_gradient,const Traits& traits);


/*!
\ingroup PkgInterpolation2Interpolation

The function `sibson_c1_interpolation` interpolates the function values and the 
gradients that are provided by functors 
following the method described in \cite s-bdnni-81. 

This function generates the interpolated function value at the point
`p` using Sibson's \f$ Z^1\f$ interpolant \cite s-bdnni-81.

If the functor `function_gradient` cannot supply the gradient of a
point, the function returns a pair where the Boolean is set to
`false`. If the interpolation was successful, the pair contains the
interpolated function value as first and `true` as second value. \pre
`norm` \f$ \neq0\f$. `function_value(p).second == true` for all points
`p` of the point/coordinate pairs in the range \f$\left[\right.\f$`first`, `beyond`\f$ \left.\right)\f$.


Parameters 
-------------- 

The template parameter `Traits` is to be 
instantiated with a model of `InterpolationTraits`. 
`ForwardIterator::value_type` is a pair associating a point to a 
(non-normalized) barycentric coordinate. `norm` is the 
normalization factor. The range \f$ \left[\right.\f$ 
`first`,`beyond`\f$ \left.\right)\f$ contains the barycentric 
coordinates for the query point `p`. The functor 
`function_value` allows to access the value of the interpolated 
function given a point. `function_gradient` allows to access the 
function gradient given a point. 

Requirements 
-------------- 

<OL> 
<LI>`Traits` is a model of the concept 
`InterpolationTraits`. 
<LI>`ForwardIterator::value_type` is a point/coordinate pair. 
Precisely `ForwardIterator::value_type::first_type` is 
equivalent to `Traits::Point_d` and 
`ForwardIterator::value_type::second_type` is equivalent to 
`Traits::FT`. 
<LI>`Functor::argument_type` must be equivalent to 
`Traits::Point_d` and `Functor::result_type` is a pair of 
the function value type and a Boolean. The function value type must 
provide a multiplication and addition operation with the type 
`Traits::FT` as well as a constructor with argument \f$ 0\f$. 
<LI>`GradFunctor::argument_type` must be equivalent to 
`Traits::Point_d` and `Functor::result_type` is a pair of 
the function's gradient type and a Boolean. The 
function gradient type must provide a multiplication operation with 
`Traits::Vector_d`. 
<LI>A model of the functor types `Functor` (resp. 
`GradFunctor`) is provided by the struct `Data_access`. It 
must be instantiated accordingly with an associative container 
(e.g. \stl `std::map`) having the point type as `key_type` 
and the function value type (resp. function gradient type) as 
`mapped_type`. 
<LI>The number type `FT` provided by `Traits` must support 
the square root operation `sqrt()`. 
</OL> 

\sa `InterpolationTraits` 
\sa `GradientFittingTraits` 
\sa CGAL::Data_access<Map> 
\sa CGAL::sibson_gradient_fitting 
\sa CGAL::linear_interpolation 
\sa CGAL::Interpolation_traits_2<K> 
\sa CGAL::Interpolation_gradient_fitting_traits_2<K> 
\sa CGAL::natural_neighbor_coordinates_2 
\sa CGAL::regular_neighbor_coordinates_2 
\sa CGAL::surface_neighbor_coordinates_3 
*/
template < class ForwardIterator, class Functor, class
GradFunctor, class Traits> std::pair< typename Functor::result_type,
bool> sibson_c1_interpolation(ForwardIterator first, ForwardIterator
beyond, const typename
std::iterator_traits<ForwardIterator>::value_type::second_type&
norm, const typename
std::iterator_traits<ForwardIterator>::value_type:: first_type& p,
Functor function_value, GradFunctor function_gradient,const Traits&
traits);

/*!
\ingroup PkgInterpolation2Interpolation

The same as `CGAL::sibson_interpolation` except that no square root
operation is needed for FT.
*/
template < class ForwardIterator, class Functor, class
GradFunctor, class Traits> typename Functor::result_type
sibson_c1_interpolation_square(ForwardIterator first,
ForwardIterator beyond, const typename
std::iterator_traits<ForwardIterator>::value_type::second_type&
norm, const typename
std::iterator_traits<ForwardIterator>::value_type::first_type& p,
Functor function_value, GradFunctor function_gradient,const
Traits& traits);

} /* namespace CGAL */

