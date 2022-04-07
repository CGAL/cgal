
namespace CGAL {

/*!
\ingroup PkgInterpolation2Ref

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
struct Data_access
  : public CGAL::cpp98::unary_function<typename Map::key_type,
                                       std::pair<typename Map::mapped_type, bool> >
{
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
If there is an entry for `p` in the container `map`, then the
pair of `map.find(p)` and `true` is returned. Otherwise, the
Boolean value of the pair is `false`.
*/
std::pair<Data_type, bool> operator()(const Key_type& p);

/// @}

}; /* end Data_access */

/*!
\ingroup PkgInterpolation2Interpolation

The function `linear_interpolation()` computes the weighted sum of the function
values which must be provided via a functor.

\tparam CoordinateInputIterator must be a model of `ForwardIterator` and
must have as value type a pair associating an entity, e.g. the `Vertex_handle` or `Point`
types of a triangulation, to a (non-normalized) barycentric coordinate.
\tparam ValueFunctor must be a functor where `ValueFunctor::argument_type` must be equivalent to
`std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` and
`ValueFunctor::result_type` is a pair of the function value type and a Boolean.
The function value type `VT` must provide an addition operator, and a multiplication operator with the type
`Traits::FT`.

A model of the functor `ValueFunctor` is provided by the struct `CGAL::Data_access` instantiated
with an associative container (e.g. `std::map`) and having:
- `std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` (the entity type) as `key_type`
- `std::iterator_traits<CoordinateInputIterator>::%value_type::second_type` (the coordinate type) as `mapped_type`.

The two template parameters must satisfy the following conditions:
- `std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` (the entity type) is equivalent to a
`ValueFunctor::argument_type`.
- `std::iterator_traits<CoordinateInputIterator>::%value_type::second_type` (the coordinate type) is a field number type
that is equivalent to `ValueFunctor::result_type::first_type`.

For example, if `CoordinateInputIterator` is an iterator with value type
`std::pair<Vertex_handle, double>`, then the `ValueFunctor` must have argument
type `Vertex_handle` (or convertible to) and return type `std::pair<double, bool>`.

\param first, beyond is the iterator range for the coordinates.
\param norm is the normalization factor.
\param value_function is a functor of type `ValueFunctor` that allows to access a pair of a
function value and a Boolean at a given entity. The Boolean indicates whether the
function value could be retrieved correctly. This function generates
the interpolated function value as the weighted sum of the values
corresponding to each entry of the entity/coordinate pairs in the range `[first, beyond)`.

\pre `norm` \f$ \neq0\f$.
\pre `first != beyond`.
\pre `value_function(p.first).second == true` for all pairs `p` in the range `[first, beyond)`.

\sa `CGAL::quadratic_interpolation()`
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa `PkgInterpolationSurfaceNeighborCoordinates3`
*/
template < class CoordinateInputIterator, class ValueFunctor >
typename ValueFunctor::result_type::first_type
linear_interpolation(CoordinateInputIterator first, CoordinateInputIterator beyond,
                     const typename std::iterator_traits<CoordinateInputIterator>::value_type::second_type& norm,
                     ValueFunctor value_function);

/*!
\ingroup PkgInterpolation2Interpolation

The function `quadratic_interpolation()` generates the
interpolated function value as the weighted sum of the values plus a
linear term in the gradient for each entity of the entity/coordinate
pairs in the range `[first, beyond)`.

\return If the interpolation was successful, the pair contains the
interpolated function value as first and `true` as second value.
Otherwise, the second value will be `false`.

\tparam Traits must be a model of `InterpolationTraits`.
Note that, contrary to some other interpolation methods, the number type `FT` provided
by `Traits` does not need to provide the square root operation.
\tparam CoordinateInputIterator must be a model of `ForwardIterator` and must have as
value type a pair associating an entity to a (non-normalized) barycentric coordinate.
More precisely, `std::iterator_traits<CoordinateInputIterator>::%value_type::first_type`
can be of the following types:
<ul>
  <li> a type equivalent to `Traits::Point_d` or `Traits::Weighted_point_d` </li>
  <li> an iterator type providing a `point()` function returning a type equivalent to `Traits::Point_d` or `Traits::Weighted_point_d`; </li>
</ul>
 and `std::iterator_traits<CoordinateInputIterator>::%value_type::second_type` must be equivalent to
`Traits::FT`.
\tparam ValueFunctor must be a functor where `ValueFunctor::argument_type` must be equivalent to
`std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` and
`ValueFunctor::result_type` is a pair of the function value type and a Boolean.
The function value type must provide a multiplication and addition operation with the type
`Traits::FT` as well as a constructor with argument `0`.
\tparam GradFunctor must be a functor where `GradFunctor::argument_type` must be equivalent to
`std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` and
`Functor::result_type` is a pair of the function's gradient type and a Boolean.
The function gradient type must provide a multiplication operation with `Traits::Vector_d`.
\tparam Point must be equivalent to `Traits::Point_d` or `Traits::Weighted_point_d`.

A model of the functor types `ValueFunctor` (resp. `GradFunctor`) is provided
by the struct `CGAL::Data_access`. It must be instantiated accordingly with an associative container
(e.g. `std::map`) having `std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` as `key_type`
and the function value type (resp. the function gradient type) as `mapped_type`.

\param first, beyond is the iterator range of the barycentric coordinates for the query point `p`.
\param norm is the normalization factor.
\param p is the point at which the interpolated function value is computed.
\param value_function is a functor that allows to access values of the interpolated function.
\param gradient_function is a functor that allows to access the function gradients.
If the functor `gradient_function` cannot supply the gradient of a point,
the function returns a pair where the Boolean is set to `false`.
\param traits is an instance of the traits class.

\pre `norm` \f$ \neq0\f$.
\pre `first != beyond`.
\pre `value_function(p.first).second == true` for pairs `p` in the range `[first, beyond)`

\sa `CGAL::linear_interpolation()`
\sa `CGAL::Interpolation_traits_2<K>`
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa `PkgInterpolationSurfaceNeighborCoordinates3`
*/
template < class CoordinateInputIterator, class ValueFunctor, class GradFunctor, class Traits, class Point >
typename ValueFunctor::result_type
quadratic_interpolation(CoordinateInputIterator first, CoordinateInputIterator beyond,
                        const typename std::iterator_traits<CoordinateInputIterator>::value_type::second_type& norm,
                        const Point& p,
                        ValueFunctor value_function,
                        GradFunctor gradient_function,
                        const Traits& traits);

/*!
\ingroup PkgInterpolation2Interpolation

The function `sibson_c1_interpolation()` generates the interpolated
function value at the point `p`, using functors for the function values
and the gradients, by applying Sibson's \f$ Z^1\f$ interpolant.

\return If the interpolation was successful, the pair contains the
interpolated function value as first and `true` as second value.
Otherwise, `false` is returned as second value.

\tparam Traits must be a model of `InterpolationTraits`.
The number type `FT` provided by `Traits` must support the square root operation `sqrt()`.
\tparam CoordinateInputIterator must be a model of `ForwardIterator` and must have as
value type a pair associating an entity to a (non-normalized) barycentric coordinate.
More precisely, `std::iterator_traits<CoordinateInputIterator>::%value_type::first_type`
can be of the following types:
<ul>
  <li> a type equivalent to `Traits::Point_d` or `Traits::Weighted_point_d` </li>
  <li> an iterator type providing a `point()` function returning a type equivalent to `Traits::Point_d` or `Traits::Weighted_point_d`; </li>
</ul>
 and `std::iterator_traits<CoordinateInputIterator>::%value_type::second_type` must be equivalent to
`Traits::FT`.
\tparam ValueFunctor must be a functor where `ValueFunctor::argument_type` must be equivalent to
`std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` and
`ValueFunctor::result_type` is a pair of the function value type and a Boolean.
The function value type must provide a multiplication and addition operation with the type
`Traits::FT` as well as a constructor with argument `0`.
\tparam GradFunctor must be a functor where `GradFunctor::argument_type` must be equivalent to
`std::iterator_traits<CoordinatetCoordinateInputIteratorIterator>::%value_type::first_type` and
`Functor::result_type` is a pair of the function's gradient type and a Boolean.
The function gradient type must provide a multiplication operation with `Traits::Vector_d`.
\tparam Point must be equivalent to `Traits::Point_d` or `Traits::Weighted_point_d`.

A model of the functor types `ValueFunctor` (resp. `GradFunctor`) is provided
by the struct `CGAL::Data_access`. It must be instantiated accordingly with an associative container
(e.g. `std::map`) having `std::iterator_traits<CoordinateInputIterator>::%value_type::first_type` as `key_type`
and the function value type (resp. the function gradient type) as `mapped_type`.

\param first, beyond is the iterator range of the barycentric coordinates for the query point `p`.
\param norm is the normalization factor.
\param p is the point at which the interpolated function value is computed.
\param value_function is a functor that allows to access values of the interpolated function.
\param gradient_function is a functor that allows to access the function gradients.
If the functor `gradient_function` cannot supply the gradient of a point,
the function returns a pair where the Boolean is set to `false`.
\param traits is an instance of the traits class.

\pre `norm` \f$ \neq0\f$.
\pre `first != beyond`.
\pre `value_function(q).second == true` for all points `q` of the point/coordinate pairs in the range `[first, beyond)`

\sa `CGAL::Interpolation_traits_2<K>`
\sa `CGAL::sibson_c1_interpolation_square()`
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3
*/
template < class CoordinateInputIterator, class ValueFunctor, class GradFunctor, class Traits, class Point >
std::pair<typename ValueFunctor::result_type, bool>
sibson_c1_interpolation(CoordinateInputIterator first, CoordinateInputIterator beyond,
                        const typename std::iterator_traits<CoordinateInputIterator>::value_type::second_type& norm,
                        const Point& p,
                        ValueFunctor value_function,
                        GradFunctor gradient_function,
                        const Traits& traits);

/*!
\ingroup PkgInterpolation2Interpolation

Same as `sibson_c1_interpolation()`, except that no square root operation is required
for the number type `Traits::FT`.

\sa `CGAL::sibson_c1_interpolation()`
*/
template < class CoordinateInputIterator, class ValueFunctor, class GradFunctor, class Traits, class Point >
typename ValueFunctor::result_type
sibson_c1_interpolation_square(CoordinateInputIterator first, CoordinateInputIterator beyond,
                               const typename std::iterator_traits<CoordinateInputIterator>::value_type::second_type& norm,
                               const Point& p,
                               ValueFunctor value_function,
                               GradFunctor gradient_function,
                               const Traits& traits);

/*!
\ingroup PkgInterpolation2Interpolation

Generates the interpolated function value computed by Farin's interpolant.

\cgalHeading{Requirements}

Same requirements as the function `sibson_c1_interpolation()`, but the
input iterator must provide random access (be a model of `RandomAccessIterator`)
and `Traits::FT` does not need to provide the square root operation.

\pre The range `[first, beyond)` contains either one or more than three elements.
The function `farin_c1_interpolation()` interpolates the function values and the
gradients that are provided by functors using the method described in \cgalCite{f-sodt-90}.

\sa `CGAL::sibson_c1_interpolation()`
*/
template < class CoordinateInputIterator, class ValueFunctor, class GradFunctor, class Traits, class Point>
typename ValueFunctor::result_type
farin_c1_interpolation(CoordinateInputIterator first, CoordinateInputIterator beyond,
                       const typename std::iterator_traits<CoordinateInputIterator>::value_type::second_type& norm,
                       const Point& p,
                       ValueFunctor value_function,
                       GradFunctor gradient_function,
                       const Traits& traits);

} /* namespace CGAL */

