namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_2` sorts iterator ranges of 
`Traits::Point_2` along a Hilbert curve by recursively subdividing 
at the median 
or the middle depending on the `PolicyTag`. 


\tparam Traits must be a model of the concept `SpatialSortingTraits_2`. 

\tparam PolicyTag is used to specify the strategy policy.
Possible values are \link CGAL::Hilbert_sort_median_policy `Hilbert_sort_median_policy` \endlink
(the default policy) or \link CGAL::Hilbert_sort_middle_policy `Hilbert_sort_middle_policy` \endlink.

\tparam ConcurrencyTag must be `Sequential_tag`, `Parallel_tag`, or `Parallel_if_available_tag. With parallelism
and TBB enabled, for the median policy up to four threads are used in parallel. 
*/
  template< typename Traits, typename PolicyTag, typename ConcurrencyTag = Sequential_tag >
class Hilbert_sort_2 {
public:

/// \name Creation 
/// @{

/*!
constructs an instance with `traits` as traits class instance. 
*/ 
Hilbert_sort_2(const Traits &traits = Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
It sorts the range `[begin, end)`. 
\tparam InputPointIterator must be a model of `RandomAccessIterator` with value type `Traits::Point_2`.
*/ 
template <class InputPointIterator>
void operator() (InputPointIterator begin, InputPointIterator end) const;

/// @}

}; /* end Hilbert_sort_2 */
} /* end namespace CGAL */
