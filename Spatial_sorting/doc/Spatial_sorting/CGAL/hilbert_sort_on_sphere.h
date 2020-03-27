namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `hilbert_sort_on_sphere()` sorts an iterator range of points that are supposed to be close to a given sphere
along a Hilbert curve on that same sphere.
If input points are not close to the input sphere, this function still works, but it might not be a good sorting function.

It sorts the range `[begin, end)` in place.

\tparam InputPointIterator must be a model of `RandomAccessIterator` and
`std::iterator_traits<InputPointIterator>::%value_type` must be convertible to `Traits::Point_3`.

\tparam Traits must be a model for concept `SpatialSortingTraits_3`.
The default traits class `Default_traits` is the kernel in which the type
`std::iterator_traits<InputPointIterator>::%value_type` is defined.

\tparam PolicyTag is used to specify the strategy policy.
Possible values are \link CGAL::Hilbert_sort_median_policy `Hilbert_sort_median_policy` \endlink
(the default policy) or \link CGAL::Hilbert_sort_middle_policy `Hilbert_sort_middle_policy` \endlink.

The input sphere is given by a squared radius and a center, parameter `sqr_radius` and parameter `center` respectively.
The default squared radius of the sphere is 1.0.
The default center of the sphere is the origin (0,0,0).

\pre `sqr_radius` greater than 0.

\cgalHeading{Implementation}

Creates an instance of `Hilbert_sort_on_sphere<Traits, PolicyTag>`,
and calls its `operator()`.

*/
template <class InputPointIterator, class Traits, class PolicyTag>
void
hilbert_sort_on_sphere( InputPointIterator begin,
                        InputPointIterator end,
                        const Traits& traits = Default_traits,
                        PolicyTag policy = Default_policy,
                        double sqr_radius = 1.0,
                        const Traits::Point_3 &center = Default_center);

} /* namespace CGAL */

