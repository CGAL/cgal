namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `spatial_sort_on_sphere()` sorts an iterator range of points in a way that
improves space locality with respect to the intrinsic metric on the sphere given as input.
Two points close in the order will be close
on the sphere, and two points close on the sphere will have a high probability
of being close in the order. The input points are supposed to be close to the input sphere.
If input points are not close to the input sphere, this function still works, but it might not be a good sorting function.

It sorts the range `[begin, end)` in place.

The input sphere is given by a squared radius and a center, parameter `sqr_radius` and parameter `center` respectively.
The default squared radius of the sphere is 1.0.
The default center of the sphere is the origin (0,0,0).

\tparam InputPointIterator must be a model of `RandomAccessIterator` and
`std::iterator_traits<InputPointIterator>::%value_type` must be convertible to
`Traits::Point_3`.

\tparam Traits must be a model for concept `SpatialSortingTraits_3`.
The default traits class `Default_traits` is the kernel in which the type
`std::iterator_traits<InputPointIterator>::%value_type` is defined.

\tparam PolicyTag is used to specify the strategy policy.
Possible values are \link CGAL::Hilbert_sort_median_policy `Hilbert_sort_median_policy` \endlink
(the default policy) or \link CGAL::Hilbert_sort_middle_policy `Hilbert_sort_middle_policy` \endlink.

\pre `sqr_radius` greater than 0.

\cgalHeading{Implementation}

Creates an instance of `Multiscale_sort<Hilbert_sort_on_sphere_3>`
where `Hilbert_sort_on_sphere_3` is an Hilbert sorting on the sphere object,
and calls its `operator()`.

The `threshold_hilbert` is the minimal size of a point set to be
subdivided recursively during Hilbert sorting, otherwise random order is used.
The `threshold_multiscale` value is the minimal size for a sample to
call the `Hilbert_sort_on_sphere_3` functor, otherwise random order is used.
The `ratio` value is used to split the original set in two subsets,
spatial sort is applied on the first subset of size
`ratio`
times the original size of the set, `Hilbert_sort_on_sphere_3` functor is applied on the
second subset.

*/
template <class InputPointIterator, class Traits, class PolicyTag>
void
spatial_sort_on_sphere( InputPointIterator begin,
                        InputPointIterator end,
                        const Traits& traits = Default_traits,
                        PolicyTag policy = Default_policy,
                        double sqr_radius = 1.0,
                        const Traits::Point_3& center = Default_center,
                        std::ptrdiff_t threshold_hilbert=default,
                        std::ptrdiff_t threshold_multiscale=default,
                        double ratio=default);

} /* namespace CGAL */

