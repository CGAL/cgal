namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `spatial_sort()` sorts an iterator range of points in a way that
improves space locality. Two points close in the order will be close
geometrically, and two points close geometrically will have a high probability
of being close in the order.

It sorts the range `[begin, end)` in place.

\tparam ConcurrencyTag enables sequential versus parallel algorithm.
Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
With parallelism enabled, sorting will be performed using up to four threads in 2D,
and up to eight threads in 3D.
Parallel sorting is available only when the median strategy policy (the default policy) is used.

\tparam InputPointIterator must be a model of `RandomAccessIterator` and
`std::iterator_traits<InputPointIterator>::%value_type` must be convertible to
`Traits::Point_2`, `Traits::Point_3`, or `Traits::Point_d`.

\tparam Traits  must be a model for concept `SpatialSortingTraits_2`,
`SpatialSortingTraits_3`, or `SpatialSortingTraits_d`.
The default traits class `Default_traits` is the kernel in which the type
`std::iterator_traits<InputPointIterator>::%value_type` is defined.

\tparam PolicyTag is used to specify the strategy policy.
Possible values are \link CGAL::Hilbert_sort_median_policy `Hilbert_sort_median_policy` \endlink
(the default policy) or \link CGAL::Hilbert_sort_middle_policy `Hilbert_sort_middle_policy` \endlink.

The default values for the thresholds and the ratio depend on the dimension.

\cgalHeading{Implementation}

Creates an instance of `Multiscale_sort<Hilbert_sort>`
where `Hilbert_sort` is an Hilbert sorting object,
and calls its `operator()`.

The `threshold_hilbert` is the minimal size of a point set to be
subdivided recursively during Hilbert sorting, otherwise random order is used.
The `threshold_multiscale` value is the minimal size for a sample to
call Hilbert sort, otherwise random order is used.
The `ratio` value is used to split the original set in two subsets,
spatial sort is applied on the first subset of size
`ratio`
times the original size of the set, Hilbert sort is applied on the
second subset.

*/
template <class ConcurrencyTag = Sequential_tag, class InputPointIterator, class Traits, class PolicyTag>
void
spatial_sort( InputPointIterator begin,
              InputPointIterator end,
              const Traits& traits = Default_traits,
              PolicyTag policy = Default_policy,
              std::ptrdiff_t threshold_hilbert=default,
              std::ptrdiff_t threshold_multiscale=default,
              double ratio=default);

} /* namespace CGAL */

