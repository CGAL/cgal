namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `hilbert_sort_on_sphere()` sorts an iterator range of points that are supposed to be close to a unit sphere
along a Hilbert curve on the unit sphere. Actually, it approximates a Hilbert curve on
the sphere by a Hilbert curve on a cube.

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

Creates an instance of `Hilbert_sort_on_sphere<Traits, PolicyTag>`,
and calls its `operator()`. 

*/
template <class RandomAccessIterator, class Traits, class PolicyTag>
void
hilbert_sort_on_sphere( RandomAccessIterator begin,
RandomAccessIterator end,
const Traits& traits = Default_traits,
PolicyTag policy = Default_policy);

} /* namespace CGAL */

