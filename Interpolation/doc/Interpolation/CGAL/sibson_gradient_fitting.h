namespace CGAL {

/*!
\defgroup PkgInterpolationSibsonGradientFitting Sibson Gradient Fitting Functions
\ingroup PkgInterpolation2Interpolation

The first function implements Sibson's gradient estimation method
based on natural neighbor coordinates \cgalCite{s-bdnni-81}.

Further functions are provided to fit the gradient for all data points
that lie inside the convex hull of the data points. One function exists
for each type of neighbor coordinates (natural and regular).

\cgalHeading{Output Format}

The output format of the functions `CGAL::sibson_gradient_fitting_nn_2()`
and `CGAL::sibson_gradient_fitting_rn_2()` can be customized using
the functor parameter of type `OutputFunctor`: this functor must have
argument type `std::pair<Tr::Vertex_handle, Traits::Vector_d>` (where `Tr`
and `Traits` are the types of the triangulation and traits passed in arguments)
but its result type is set as desired by the user.

See also \link PkgInterpolationNaturalNeighborCoordinates2 natural neighbor coordinates functions\endlink,
which use the same mechanism to allow flexible output.

\sa `CGAL::Interpolation_gradient_fitting_traits_2<K>`
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa `PkgInterpolationSurfaceNeighborCoordinates3`
\sa `PkgInterpolation2Interpolation`
*/
/// @{

/*!
Estimates the gradient of a function at a query point.

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
\tparam Traits must be a model of the concept `GradientFittingTraits`.
\tparam Point must be equivalent to `Traits::Point_d` or `Traits::Weighted_point_d`.

An example of compatible types for `CoordinateInputIterator` and `ValueFunctor`
is a forward iterator with value type `std::pair<Vertex_handle, Traits::FT>`
and a functor with argument type `Vertex_handle`.

\param first, beyond is the iterator range of the neighbor coordinates for the query point `p`.
\param norm is the normalization factor.
\param p is the query point
\param value_function is a functor that allows to access values of the interpolated function.
\param fn is the value of the function at `p`.
\param traits is an instance of the gradient fitting traits class.

\return the gradient at the point `p`.
*/
template < class CoordinateInputIterator, class ValueFunctor, class Traits, class Point >
typename Traits::Vector_d
sibson_gradient_fitting(CoordinateInputIterator first, CoordinateInputIterator beyond,
                        const typename std::iterator_traits<CoordinateInputIterator>::value_type::second_type& norm,
                        const Point& p,
                        ValueFunctor value_function,
                        const typename ValueFunctor::result_type::first_type fn,
                        const Traits& traits);

/*!
Estimates the function gradients at all vertices of the Delaunay triangulation `dt`
that lie inside the convex hull, using the coordinates computed by the
function \ref PkgInterpolationNaturalNeighborCoordinates2.

\tparam Dt must be of type `Delaunay_triangulation_2<Dt_Traits, Tds>`.
        `Dt_Traits` must be a model of the concepts `DelaunayTriangulationTraits_2` and `PolygonTraits_2`.
\tparam GradientOutputIterator must be a model of `OutputIterator` with value type
        `OutputFunctor::result_type`.
\tparam OutputFunctor must be a functor with argument type `std::pair<Dt::Vertex_handle, Traits::Vector_d>`.
        Note that the result type of the functor is not specified and is chosen by users to fit their needs.
\tparam ValueFunctor must be a functor where:
- `ValueFunctor::argument_type` must be either `Dt::Vertex_handle` or `Dt::Point`.
- `ValueFunctor::result_type` is a pair of the function value type and a Boolean.
The function value type must provide a multiplication and addition operation with the type
`Traits::FT` as well as a constructor with argument `0`.
\tparam Traits must be a model of `GradientFittingTraits`.

\param dt is the Delaunay triangulation.
\param out is an object of type `GradientOutputIterator`.
\param fct is an object of type `OutputFunctor`.
\param value_function is a functor of type `ValueFunctor` that gives access to
the values of the function at points of the triangulation.
\param traits is an instance of the gradient fitting traits class.

\return An output iterator with value type `OutputFunctor::result_type`.
The sequence is written starting at `out`. The function returns an
iterator that is placed past-the-end of the resulting sequence.
*/
template < class Dt, class GradientOutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
GradientOutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             GradientOutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits);

/*!
Estimates the function gradients at all vertices of `rt` that lie
inside the convex hull using the coordinates computed by the
functions \ref PkgInterpolationRegularNeighborCoordinates2.

\tparam Rt must be of type `Regular_triangulation_2<Rt_Traits, Tds>`.
        `Rt_Traits` must be a model of the concepts `RegularTriangulationTraits_2` and `PolygonTraits_2`.
\tparam GradientOutputIterator must be a model of `OutputIterator` with value type
        `OutputFunctor::result_type`.
\tparam OutputFunctor must be a functor with argument type `std::pair<Rt::Vertex_handle, Traits::Vector_d>`.
        Note that the result type of the functor is not specified and is chosen by users to fit their needs.
\tparam ValueFunctor must be a functor where:
- `ValueFunctor::argument_type` must be either `Rt::Vertex_handle` or `Rt::Weighted_point`.
- `ValueFunctor::result_type` is a pair of the function value type and a Boolean.
The function value type must provide a multiplication and addition operation with the type
`Traits::FT` as well as a constructor with argument `0`.
\tparam Traits must be a model of `GradientFittingTraits`.

\param rt is the regular triangulation.
\param out is an object of type `GradientOutputIterator`.
\param fct is an object of type `OutputFunctor`.
\param value_function is a functor of type `ValueFunctor` that gives access to
the values of the function at points of the triangulation.
\param traits is an instance of the gradient fitting traits class.

\return An output iterator with value type `OutputFunctor::result_type`.
The sequence is written starting at `out`. The function returns an
iterator that is placed past-the-end of the resulting sequence.
*/
template < class Rt, class GradientOutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
GradientOutputIterator
sibson_gradient_fitting_rn_2(const Rt& rt,
                             GradientOutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits);

/// @}

} /* namespace CGAL */

