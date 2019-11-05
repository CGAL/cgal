namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `hilbert_sort()` sorts an iterator range of points 
along a Hilbert curve. 

It sorts the range `[begin, end)` in place. 


\tparam ConcurrencyTag With `ConcurrencyTag = Parallel_tag` and with TBB found, the sorting will be 
done in 4 threads in 2D, and in 8 threads in 3D with the default policy.

\tparam RandomAccessIterator 
`std::iterator_traits<RandomAccessIterator>::%value_type` must be convertible to 
`Traits::Point_2`, `Traits::Point_3`, or `Traits::Point_d`. 

\tparam Traits must be a model for concept `SpatialSortingTraits_2`, 
`SpatialSortingTraits_3`, or `SpatialSortingTraits_d`. 
The default traits class `Default_traits` is the kernel in which the type 
`std::iterator_traits<RandomAccessIterator>::%value_type` is defined. 

\tparam PolicyTag The default policy is `Hilbert_sort_median_policy()` and the 
other option is `Hilbert_sort_middle_policy()`. 

\cgalHeading{Requirements}

<OL> 

</OL> 

\cgalHeading{Implementation}

Creates an instance of 
`Hilbert_sort_2<Traits, PolicyTag>`, 
`Hilbert_sort_3<Traits, PolicyTag>`, or 
`Hilbert_sort_d<Traits, PolicyTag>` 
and calls its `operator()`. 

*/
template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Traits, class PolicyTag>
void
hilbert_sort( RandomAccessIterator begin,
              RandomAccessIterator end,
              const Traits& traits = Default_traits,
              PolicyTag policy = Default_policy);

} /* namespace CGAL */

