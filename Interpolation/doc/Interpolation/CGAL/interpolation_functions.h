
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

\cgalHeading{Parameters}

The class 
`Data_access` has the container type `Map` as template parameter. 

*/
template< typename Map >
struct Data_access : public std::unary_function< typename Map::key_type,
		     std::pair< typename Map::mapped_type, bool> > {
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

generates the interpolated function value computed by Farin's interpolant.

\pre `norm` \f$ \neq0\f$. `function_value(p).second == true` for all points `p` of the point/coordinate pairs in the range `[first, beyond)`.
\pre The range `[first, beyond)` contains either one or more than three elements.
The function `farin_c1_interpolation()` interpolates the function values and the 
gradients that are provided by functors using the method described in \cgalCite{f-sodt-90}. 

\cgalHeading{Parameters}

The value type of `RandomAccessIterator` is a pair 
associating a point to a (non-normalized) barycentric coordinate. See 
`sibson_c1_interpolation()` for the other parameters. 

\cgalHeading{Requirements}

Same requirements as for `sibson_c1_interpolation()` only the 
iterator must provide random access and `Traits::FT` does not need 
to provide the square root operation. 

\sa `CGAL::Data_access<Map>` 
\sa `CGAL::linear_interpolation()` 
\sa `CGAL::sibson_c1_interpolation()` 
\sa PkgInterpolationSibsonGradientFitting
\sa `CGAL::Interpolation_traits_2<K>`
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3

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

The function `linear_interpolation()` computes the weighted sum of the function 
values which must be provided via a functor. 

\tparam  ForwardIterator must have as value type a pair associating a point to a 
(non-normalized) barycentric coordinate, that is
`std::iterator_traits<ForwardIterator>::%value_type::first_type` is equivalent to a 
point and `std::iterator_traits<ForwardIterator>::%value_type::second_type` is a 
field number type. 
\tparam Functor The type `Functor::argument_type` must be equivalent to 
`std::iterator_traits<ForwardIterator>::%value_type::first_type` and 
`Functor::result_type` is a pair of the function value type 
and a Boolean value. The function value type must provide a 
multiplication and addition operation with the field number type 
`std::iterator_traits<ForwardIterator>::%value_type::second_type` and a constructor 
with argument `0`. 

A model of the functor is provided by the 
struct `Data_access`. It must be instantiated accordingly with 
an associative container (e.g. `std::map`) having the 
point type as `key_type` and the function value type as 
`mapped_type`.

\param first,beyond are the iterator range for the weighted input points.
\param norm is the normalization factor. `norm` \f$ \neq0\f$.
\param function_values  is a functor that allows to access a pair of a
function value and a Boolean for a given point. The Boolean indicates whether the
function value could be retrieved correctly. This function generates
the interpolated function value as the weighted sum of the values
corresponding to each point of the point/coordinate pairs in the
range `[first, beyond)`. 
 `function_value(q).second == true` for all points `q` of the point/coordinate pairs in the range `[first, beyond)`.


\sa `CGAL::Data_access<Map>` 
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3

*/
template < class ForwardIterator, class Functor> typename
Functor::result_type::first_type linear_interpolation(ForwardIterator first,
ForwardIterator beyond, const typename
std::iterator_traits<ForwardIterator>::value_type::second_type&
norm, Functor function_values);


/*!
\ingroup PkgInterpolation2Interpolation

The function `quadratic_interpolation()` generates the
interpolated function value as the weighted sum of the values plus a
linear term in the gradient for each point of the point/coordinate
pairs in the range `[first, beyond)`. 

\cgalHeading{Parameters and Template Parameters}

The same as for `sibson_c1_interpolation()` only that `Traits::FT` does not need 
to provide the square root operation. 

\sa `InterpolationTraits` 
\sa `GradientFittingTraits` 
\sa `CGAL::Data_access<Map>` 
\sa PkgInterpolationSibsonGradientFitting
\sa `CGAL::linear_interpolation()`
\sa `CGAL::Interpolation_traits_2<K>` 
\sa `CGAL::Interpolation_gradient_fitting_traits_2<K>` 
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3
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

The function `sibson_c1_interpolation()` generates the interpolated
function value at the point `p`, using functors for the function values 
and the gradients, by applying Sibson's \f$ Z^1\f$ interpolant. 

If the functor `function_gradient` cannot supply the gradient of a
point, the function returns a pair where the Boolean is set to
`false`. If the interpolation was successful, the pair contains the
interpolated function value as first and `true` as second value. 


\tparam Traits must be a model of `InterpolationTraits`. 
The number type `FT` provided by `Traits` must support 
the square root operation `sqrt()`. 
\tparam ForwardIterator must have as value type a pair associating a point to a 
(non-normalized) barycentric coordinate.
More precisely, `std::iterator_traits<ForwardIterator>::%value_type::first_type` is 
equivalent to `Traits::Point_d` and 
`std::iterator_traits<ForwardIterator>::%value_type::second_type` is equivalent to 
`Traits::FT`. 
 
\tparam Functor must be a functor where `Functor::argument_type` must be equivalent to 
`Traits::Point_d` and `Functor::result_type` is a pair of 
the function value type and a Boolean. The function value type must 
provide a multiplication and addition operation with the type 
`Traits::FT` as well as a constructor with argument `0`.

\tparam GradFunctor must be a functor where `GradFunctor::argument_type` must be equivalent to 
`Traits::Point_d` and `Functor::result_type` is a pair of 
the function's gradient type and a Boolean. The 
function gradient type must provide a multiplication operation with 
`Traits::Vector_d`. 

A model of the functor types `Functor` (resp. 
`GradFunctor`) is provided by the struct `Data_access`. It 
must be instantiated accordingly with an associative container 
(e.g. `std::map`) having the point type as `key_type` 
and the function value type (resp. function gradient type) as 
`mapped_type`. 

\param first,beyond is the iterator range of the barycentric 
coordinates for the query point `p`. 
\param norm is the normalization factor. `norm` \f$ \neq0\f$.
\param 	p is the point at which the interpolated function value is generated
\param function_value is a functor that allows to access the value of the interpolated 
function given a point. `function_value(q).second == true` for all points
`q` of the point/coordinate pairs in the range `[first, beyond)`
\param function_gradient is a functor that allows to access the 
function gradient given a point. 
\param traits is an instance of the traits class

\sa `InterpolationTraits` 
\sa `GradientFittingTraits` 
\sa `CGAL::Data_access<Map>` 
\sa PkgInterpolationSibsonGradientFitting
\sa `CGAL::linear_interpolation()`
\sa `CGAL::Interpolation_traits_2<K>` 
\sa `CGAL::Interpolation_gradient_fitting_traits_2<K>` 
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3
*/
template < class ForwardIterator, class Functor, class
GradFunctor, class Traits> std::pair< typename Functor::result_type,
bool> sibson_c1_interpolation(ForwardIterator first, ForwardIterator
beyond, const typename
std::iterator_traits<ForwardIterator>::value_type::second_type&
norm, const typename
std::iterator_traits<ForwardIterator>::value_type::first_type& p,
Functor function_value, GradFunctor function_gradient,const Traits&
traits);

/*!
\ingroup PkgInterpolation2Interpolation

The same as `sibson_c1_interpolation()` except that no square root
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

