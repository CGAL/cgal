namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `hilbert_sort()` sorts an iterator range of points
along a Hilbert curve.

It sorts the range `[begin, end)` in place.

\tparam ConcurrencyTag enables sequential versus parallel algorithm.
Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
With parallelism enabled, sorting will be performed using up to four threads in 2D,
and up to eight threads in 3D.
Parallel sorting is available only when the median strategy policy (the default policy) is used.

\tparam InputPointIterator must be a model of `RandomAccessIterator` and
`std::iterator_traits<InputPointIterator>::%value_type` must be convertible to
`Traits::Point_2`, `Traits::Point_3`, or `Traits::Point_d`.

\tparam Traits must be a model for concept `SpatialSortingTraits_2`,
`SpatialSortingTraits_3`, or `SpatialSortingTraits_d`.
The default traits class `Default_traits` is the kernel in which the type
`std::iterator_traits<InputPointIterator>::%value_type` is defined.

\tparam PolicyTag is used to specify the strategy policy.
Possible values are \link CGAL::Hilbert_sort_median_policy `Hilbert_sort_median_policy` \endlink
(the default policy) or \link CGAL::Hilbert_sort_middle_policy `Hilbert_sort_middle_policy` \endlink.

\cgalHeading{Implementation}

Creates an instance of
`Hilbert_sort_2<Traits, PolicyTag>`,
`Hilbert_sort_3<Traits, PolicyTag>`, or
`Hilbert_sort_d<Traits, PolicyTag>`
and calls its `operator()`.

*/
template <class ConcurrencyTag = Sequential_tag, class InputPointIterator, class Traits, class PolicyTag>
void
hilbert_sort( InputPointIterator begin,
              InputPointIterator end,
              const Traits& traits = Default_traits,
              PolicyTag policy = Default_policy);

} /* namespace CGAL */

