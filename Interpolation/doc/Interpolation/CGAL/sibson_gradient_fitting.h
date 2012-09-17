namespace CGAL {

/*!
\defgroup sibson_gradient_fitting sibson_gradient_fitting
\ingroup PkgInterpolation2Interpolation

The function `sibson_gradient_fitting` approximates the gradient of a
function at a point `p` given natural neighbor coordinates for `p` and
its neighbors' function values. The approximation method is described
in \cite s-bdnni-81. Further functions are provided to fit the
gradient for all data points that lie inside the convex hull of the
data points. One function exists for each type of natural neighbor
coordinates.

### Requirements ###

<OL> 
<LI>`ForwardIterator::value_type` is a pair of point/coordinate 
value, thus `ForwardIterator::value_type::first_type` is 
equivalent to a point and 
`ForwardIterator::value_type::second_type` is a 
number type. 
<LI>`Functor::argument_type` must be equivalent to 
`ForwardIterator::value_type::first_type` and 
`Functor::result_type` is the function value type. It must 
provide a multiplication and addition operation with the type 
`ForwardIterator::value_type::second_type`. 
<LI>`Traits` is a model of the concept 
`GradientFittingTraits`. 
</OL> 

\sa CGAL::linear_interpolation 
\sa CGAL::sibson_c1_interpolation 
\sa CGAL::farin_c1_interpolation 
\sa CGAL::quadratic_interpolation 
\sa `CGAL::Interpolation_gradient_fitting_traits_2<K>`
\sa CGAL::natural_neighbor_coordinates_2 
\sa CGAL::regular_neighbor_coordinates_2 
\sa CGAL::surface_neighbor_coordinates_3 

### Implementation ###

This function implements Sibson's gradient 
estimation method based on natural neighbor coordinates 
\cite s-bdnni-81. 

*/
/// @{

/*!
This function estimates the
gradient of a function at the point `p` given natural neighbor
coordinates of `p` in the range \f$ \left[\right.\f$ `first`,
`beyond`\f$ \left.\right)\f$ and the function values of the neighbors
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
function `CGAL::natural_neighbor_coordinates_2`.
`OutputIterator::value_type` is a pair associating a point to a
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
function `CGAL::regular_neighbor_coordinates_2`.
`OutputIterator::value_type` is a pair associating a point to a
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

