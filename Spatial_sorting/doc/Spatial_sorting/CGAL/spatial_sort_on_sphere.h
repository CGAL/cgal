namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `spatial_sort_on_sphere()` sorts an iterator range of points in a way that 
improves space locality. Two points close in the order will be close 
on the sphere, and two points close on the sphere will have a high probability 
of being close in the order. 

sorts the range `[begin, end)` in place. 

The default traits class `Default_traits` is the kernel in which the type 
`std::iterator_traits<RandomAccessIterator>::%value_type` is defined. 

The default policy is `Hilbert_sort_median_policy()` and the 
other option is `Hilbert_sort_middle_policy()`. 

\cgalHeading{Requirements}

<OL> 
<LI>`std::iterator_traits<RandomAccessIterator>::%value_type` is convertible to 
`Traits::Point_3`. 
<LI>`Traits` is a model for concept `SpatialSortingTraits_3`. 
</OL> 

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
template <class RandomAccessIterator, class Traits, class PolicyTag>
void
spatial_sort_on_sphere( RandomAccessIterator begin,
RandomAccessIterator end,
const Traits& traits = Default_traits,
PolicyTag policy = Default_policy,
std::ptrdiff_t threshold_hilbert=default,
std::ptrdiff_t threshold_multiscale=default,
double ratio=default);

} /* namespace CGAL */

