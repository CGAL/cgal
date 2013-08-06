namespace CGAL {

/*!
\defgroup PkgInterpolationSibsonGradientFitting Sibson gradient fitting functions
\ingroup PkgInterpolation2Interpolation

These functions approximate the gradient of a
function at a point `p` given natural neighbor coordinates for `p` and
its neighbors' function values. The approximation method is described
in \cgalCite{s-bdnni-81}. Further functions are provided to fit the
gradient for all data points that lie inside the convex hull of the
data points. One function exists for each type of natural neighbor
coordinates.

\cgalHeading{Requirements}

<OL> 
<LI>The value type of `ForwardIterator` is a pair of point/coordinate 
value, thus `std::iterator_traits<ForwardIterator>::%value_type::first_type` is 
equivalent to a point and 
`std::iterator_traits<ForwardIterator>::%value_type::second_type` is a 
number type. 
<LI>`Functor::argument_type` must be equivalent to 
`std::iterator_traits<ForwardIterator>::%value_type::first_type` and 
`Functor::result_type` is the function value type. It must 
provide a multiplication and addition operation with the type 
`std::iterator_traits<ForwardIterator>::%value_type::second_type`. 
<LI>`Traits` is a model of the concept 
`GradientFittingTraits`. 
</OL> 

\sa `CGAL::linear_interpolation()`
\sa `CGAL::sibson_c1_interpolation()` 
\sa `CGAL::farin_c1_interpolation()`
\sa `CGAL::quadratic_interpolation()` 
\sa `CGAL::Interpolation_gradient_fitting_traits_2<K>`
\sa `PkgInterpolationNaturalNeighborCoordinates2`
\sa `PkgInterpolationRegularNeighborCoordinates2`
\sa PkgInterpolationSurfaceNeighborCoordinates3

\cgalHeading{Implementation}

This function implements Sibson's gradient 
estimation method based on natural neighbor coordinates 
\cgalCite{s-bdnni-81}. 

*/
/// @{

/*!
estimates the
gradient of a function at the point `p` given natural neighbor
coordinates of `p` in the range `[first, beyond)` and the function values of the neighbors
provided by the functor `f`. `norm` is the normalization
factor of the barycentric coordinates.
*/
template < class ForwardIterator, class Functor, class
Traits> typename Traits::Vector_d
sibson_gradient_fitting(ForwardIterator first, ForwardIterator
beyond, const typename
std::iterator_traits<ForwardIterator>::value_type::second_type&
norm, const typename
std::iterator_traits<ForwardIterator>::value_type::first_type& p,
Functor f, const Traits& traits);

/*!
estimates the function gradients at all vertices of `dt` that lie
inside the convex hull using the coordinates computed by the
function `PkgInterpolationNaturalNeighborCoordinates2`.
The value type of `OutputIterator` is a pair associating a point to a
vector. The sequence of point/gradient pairs computed by this
function is placed starting at `out`. The function returns an
iterator that is placed past-the-end of the resulting sequence. The
requirements are the same as above. The template class `Dt` must
be equivalent to `Delaunay_triangulation_2<Gt, Tds>`.
*/
template < class Dt, class OutputIterator, class Functor,
class Traits> OutputIterator sibson_gradient_fitting_nn_2(const Dt&
dt, OutputIterator out, Functor f, const Traits& traits);

/*!
estimates the function gradients at all vertices of `rt` that lie
inside the convex hull using the coordinates computed by the
functions `PkgInterpolationRegularNeighborCoordinates2`.
The value type of `OutputIterator` is a pair associating a point to a
vector. The sequence of point/gradient pairs computed by this
function is placed starting at `out`. The function returns an
iterator that is placed past-the-end of the resulting sequence. The
requirements are the same as above. The template class `Rt` must
be equivalent to `Regular_triangulation_2<Gt, Tds>`.
*/
template < class Rt, class OutputIterator, class Functor,
class Traits> OutputIterator sibson_gradient_fitting_rn_2(const Rt&
rt, OutputIterator out, Functor f, const Traits& traits);

/// @}

} /* namespace CGAL */

